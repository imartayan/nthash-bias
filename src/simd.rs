use seq_hash::packed_seq::u32x8 as S;

/// Count leading zeros for each u32 lane of a `u32x8`.
#[inline(always)]
pub fn lzcnt_u32x8(v: S) -> S {
    #[cfg(target_feature = "avx2")]
    {
        lzcnt_avx2(v)
    }
    #[cfg(target_feature = "neon")]
    {
        lzcnt_neon(v)
    }
    #[cfg(not(any(target_feature = "avx2", target_feature = "neon")))]
    {
        let a = v.to_array();
        S::from([
            a[0].leading_zeros(),
            a[1].leading_zeros(),
            a[2].leading_zeros(),
            a[3].leading_zeros(),
            a[4].leading_zeros(),
            a[5].leading_zeros(),
            a[6].leading_zeros(),
            a[7].leading_zeros(),
        ])
    }
}

/// AVX2: binary search for the position of the highest set bit.
///
/// Each round checks whether the upper half of the remaining range is zero.
/// If so, it adds the half-width to the count and shifts x left to examine
/// the lower half next. After 5 rounds (16 → 8 → 4 → 2 → 1) we have lzcnt
/// for all non-zero lanes. The x == 0 case would produce 31 from the rounds
/// alone, so a final fixup adds 1 where the original value was zero.
#[cfg(target_feature = "avx2")]
#[inline(always)]
fn lzcnt_avx2(v: S) -> S {
    use std::arch::x86_64::*;
    use std::mem::transmute;
    unsafe {
        let orig: __m256i = transmute(v);
        let mut x = orig;
        let zero = _mm256_setzero_si256();
        let mut n = zero;

        macro_rules! round {
            ($shift:literal, $bits:literal) => {{
                // mask = all-ones in lanes where the upper `$shift` bits are zero
                let mask = _mm256_cmpeq_epi32(_mm256_srli_epi32(x, $shift), zero);
                n = _mm256_add_epi32(n, _mm256_and_si256(mask, _mm256_set1_epi32($bits)));
                x = _mm256_blendv_epi8(x, _mm256_slli_epi32(x, $bits), mask);
            }};
        }

        round!(16, 16);
        round!(24, 8);
        round!(28, 4);
        round!(30, 2);
        round!(31, 1);

        // Lanes that were zero yield n=31 above; bump them to 32.
        let zero_mask = _mm256_cmpeq_epi32(orig, zero);
        n = _mm256_add_epi32(n, _mm256_and_si256(zero_mask, _mm256_set1_epi32(1)));

        transmute(n)
    }
}

#[cfg(target_feature = "neon")]
#[inline(always)]
fn lzcnt_neon(v: S) -> S {
    use std::arch::aarch64::{uint32x4_t, vclzq_u32};
    use std::mem::transmute;
    unsafe {
        let [a, b]: [uint32x4_t; 2] = transmute(v);
        transmute([vclzq_u32(a), vclzq_u32(b)])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lzcnt_u32x8() {
        // Values chosen to cover: 0, powers of 2, all-ones, near-boundary cases.
        let input = S::from([
            0u32,
            1,
            2,
            0x0000_8000,
            0x4000_0000,
            0x7FFF_FFFF,
            0x8000_0000,
            0xFFFF_FFFF,
        ]);
        let expected = S::from([32u32, 31, 30, 16, 1, 1, 0, 0]);
        assert_eq!(lzcnt_u32x8(input), expected);
    }
}
