//! A small, fast, non-cryptographic string hash (FNV-1a), used wherever a
//! command needs to compare or deduplicate read IDs across possibly
//! hundreds of millions of records without keeping every ID string in
//! memory -- `fastx deinterleave` (layout detection) and `fastx cat`
//! (duplicate-ID detection) both hash `FastqRecord::base_id()` down to a
//! single `u64` instead of storing/comparing full strings.
//!
//! FNV-1a is not collision-proof, but for read IDs (structured, mostly
//! distinct strings) the odds of an accidental 64-bit collision are
//! astronomically small next to the actual failure modes these checks are
//! built to catch (a whole file duplicated, or a merged file that isn't
//! interleaved at all).

use std::hash::{BuildHasherDefault, Hasher};

const FNV_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01B3;

pub fn fnv1a(bytes: &[u8]) -> u64 {
    fnv1a_from(FNV_OFFSET_BASIS, bytes)
}

/// FNV-1a continued from an existing state, so a hash can be built from
/// several pieces.
fn fnv1a_from(mut hash: u64, bytes: &[u8]) -> u64 {
    for &b in bytes {
        hash ^= b as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

/// A hasher for keys that are *already* hashes.
///
/// `fastx cat` keys its seen-set on `fnv1a` values of read IDs. The standard
/// library's default hasher would run SipHash over each of those u64s --
/// hashing a hash, once per read, on inputs with hundreds of millions of
/// them. This takes the u64 as given and applies only splitmix64's finalizer:
/// a few multiplies and shifts instead of a full SipHash.
///
/// The finalizer is not decoration. A hash table picks its bucket from the
/// low bits of the hash and its control byte from the top few; FNV-1a's low
/// bits are its least-mixed, since each input byte reaches them through a
/// single multiply. Feeding raw FNV values to the table would cluster
/// buckets and undo the speedup. Mixing first makes both ends uniform.
///
/// Keys that are not single u64s still hash correctly -- `write` folds bytes
/// in with FNV -- they just gain nothing from this hasher.
#[derive(Default)]
pub struct IdHasher(u64);

impl Hasher for IdHasher {
    fn write_u64(&mut self, n: u64) {
        self.0 = n;
    }

    fn write(&mut self, bytes: &[u8]) {
        self.0 = fnv1a_from(self.0, bytes);
    }

    fn finish(&self) -> u64 {
        // splitmix64's finalizer; the same mixing `common::rng` uses to turn
        // a counter into a well-distributed draw.
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// Use as `HashSet<u64, BuildIdHasher>`.
pub type BuildIdHasher = BuildHasherDefault<IdHasher>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_and_sensitive_to_input() {
        assert_eq!(fnv1a(b"READ_1"), fnv1a(b"READ_1"));
        assert_ne!(fnv1a(b"READ_1"), fnv1a(b"READ_2"));
        assert_ne!(fnv1a(b""), fnv1a(b"a"));
    }

    #[test]
    fn id_hasher_passes_u64_keys_through_a_finalizer() {
        use std::hash::BuildHasher;
        let build = BuildIdHasher::default();
        let of = |n: u64| build.hash_one(n);

        assert_eq!(of(42), of(42), "must be deterministic");
        assert_ne!(of(42), of(43));
        // Raw pass-through would leave the table clustering on FNV's weakest
        // bits, so the finalizer has to actually change the value.
        assert_ne!(of(42), 42);
    }

    #[test]
    fn id_hasher_spreads_sequential_and_fnv_keys_across_buckets() {
        use std::hash::BuildHasher;
        let build = BuildIdHasher::default();
        // Low 12 bits decide the bucket for a table of a few thousand slots.
        let buckets = |keys: Vec<u64>| {
            keys.into_iter()
                .map(|k| build.hash_one(k) & 0xfff)
                .collect::<std::collections::HashSet<_>>()
                .len()
        };
        let sequential = buckets((0..4096).collect());
        let fnv = buckets(
            (0..4096)
                .map(|i| fnv1a(format!("READ_{i:08}").as_bytes()))
                .collect(),
        );
        // ~4096 keys into 4096 buckets fills about 63% of them if the spread
        // is uniform; anything clustered lands far below that.
        assert!(sequential > 2200, "sequential keys clustered: {sequential}");
        assert!(fnv > 2200, "FNV keys clustered: {fnv}");
    }

    #[test]
    fn used_as_a_hash_set_hasher_it_still_behaves_like_a_set() {
        let mut seen: std::collections::HashSet<u64, BuildIdHasher> = Default::default();
        for i in 0..10_000u64 {
            assert!(seen.insert(fnv1a(&i.to_le_bytes())), "first insert is new");
        }
        for i in 0..10_000u64 {
            assert!(
                !seen.insert(fnv1a(&i.to_le_bytes())),
                "second insert is a repeat"
            );
        }
        assert_eq!(seen.len(), 10_000);
    }

    #[test]
    fn no_collisions_in_a_small_structured_set() {
        use std::collections::HashSet;
        let hashes: HashSet<u64> = (0..100_000)
            .map(|i| fnv1a(format!("READ_{i:08}").as_bytes()))
            .collect();
        assert_eq!(hashes.len(), 100_000);
    }
}
