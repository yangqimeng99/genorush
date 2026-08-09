//! Small dependency-free SplitMix64-based RNG, used two ways by `fastx sample`:
//! - stateless `deterministic_f64(seed, index)`: a pure function of (seed, index),
//!   so any record's keep/discard draw can be computed independently of every
//!   other record — the basis for parallel, order-independent proportion sampling.
//! - stateful `SplitMix64`: a conventional sequential PRNG, used by the reservoir
//!   sampler where each draw depends on the previous stream position.

const GOLDEN_GAMMA: u64 = 0x9E37_79B9_7F4A_7C15;

fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Deterministic, order-independent pseudorandom u64 for a given (seed, index) pair.
pub fn deterministic_u64(seed: u64, index: u64) -> u64 {
    mix64(seed.wrapping_add(index.wrapping_mul(GOLDEN_GAMMA)))
}

/// Deterministic, order-independent uniform f64 in [0, 1) for a given (seed, index) pair.
pub fn deterministic_f64(seed: u64, index: u64) -> f64 {
    (deterministic_u64(seed, index) >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(GOLDEN_GAMMA);
        mix64(self.state)
    }

    /// Unbiased uniform integer in [0, bound) via Lemire's method.
    pub fn next_below(&mut self, bound: u64) -> u64 {
        debug_assert!(bound > 0);
        let mut x = self.next_u64();
        let mut m = (x as u128) * (bound as u128);
        let mut l = m as u64;
        if l < bound {
            let t = bound.wrapping_neg() % bound;
            while l < t {
                x = self.next_u64();
                m = (x as u128) * (bound as u128);
                l = m as u64;
            }
        }
        (m >> 64) as u64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_f64_is_pure_and_in_range() {
        for seed in [0u64, 1, 42, u64::MAX] {
            for index in [0u64, 1, 2, 1_000_000, u64::MAX] {
                let a = deterministic_f64(seed, index);
                let b = deterministic_f64(seed, index);
                assert_eq!(a, b, "must be a pure function of (seed, index)");
                assert!((0.0..1.0).contains(&a), "out of range: {a}");
            }
        }
    }

    #[test]
    fn deterministic_f64_differs_across_index_and_seed() {
        let a = deterministic_f64(42, 0);
        let b = deterministic_f64(42, 1);
        let c = deterministic_f64(7, 0);
        assert_ne!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn deterministic_f64_mean_is_roughly_uniform() {
        let seed = 123;
        let n = 200_000u64;
        let sum: f64 = (0..n).map(|i| deterministic_f64(seed, i)).sum();
        let mean = sum / n as f64;
        assert!((0.48..0.52).contains(&mean), "sample mean {mean} looks biased");
    }

    #[test]
    fn next_below_is_in_bounds_and_covers_range() {
        let mut rng = SplitMix64::new(99);
        let bound = 7u64;
        let mut seen = std::collections::HashSet::new();
        for _ in 0..10_000 {
            let x = rng.next_below(bound);
            assert!(x < bound);
            seen.insert(x);
        }
        assert_eq!(seen.len() as u64, bound, "did not observe all values in [0, bound)");
    }

    #[test]
    fn same_seed_same_stream() {
        let mut a = SplitMix64::new(2024);
        let mut b = SplitMix64::new(2024);
        for _ in 0..1000 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }
}
