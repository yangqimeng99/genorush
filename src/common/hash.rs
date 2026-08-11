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

const FNV_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01B3;

pub fn fnv1a(bytes: &[u8]) -> u64 {
    let mut hash = FNV_OFFSET_BASIS;
    for &b in bytes {
        hash ^= b as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

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
    fn no_collisions_in_a_small_structured_set() {
        use std::collections::HashSet;
        let hashes: HashSet<u64> = (0..100_000)
            .map(|i| fnv1a(format!("READ_{i:08}").as_bytes()))
            .collect();
        assert_eq!(hashes.len(), 100_000);
    }
}
