use seq_hash::packed_seq::{AsciiSeqVec, PackedSeqVec, Seq, SeqVec, u32x8};
use seq_hash::{KmerHasher, MulHasher, NtHasher};

use core::hash::BuildHasher;
use core::hint::black_box;
use std::time::Instant;

const LEN: usize = 100_000_000;
const K: usize = 31;
const MASK: u64 = (1 << (2 * K)) - 1;

fn main() {
    let ascii: Vec<u8> = AsciiSeqVec::random(LEN).seq;
    let packed = PackedSeqVec::from_ascii(&ascii);

    // ascii iter + rapidhash
    {
        let t = Instant::now();
        let checksum = ascii.windows(K).fold(0u64, |acc, kmer| {
            acc.wrapping_add(rapidhash::fast::GlobalState::default().hash_one(kmer))
        });
        let elapsed = t.elapsed().as_secs_f64();
        black_box(checksum);
        eprintln!(
            "ascii seq + rapidhash:\t{:.2} GB/s",
            LEN as f64 / 1e9 / elapsed
        );
    }

    // packed-seq iter + fxhash
    {
        let t = Instant::now();
        let mut kmer: u64 = 0;
        packed.as_slice().iter_bp().take(K - 1).for_each(|bp| {
            kmer = (kmer << 2) | bp as u64;
        });
        let checksum = packed
            .as_slice()
            .iter_bp()
            .skip(K - 1)
            .fold(0u64, |acc, bp| {
                kmer = ((kmer << 2) | bp as u64) & MASK;
                acc.wrapping_add(rustc_hash::FxBuildHasher.hash_one(kmer))
            });
        let elapsed = t.elapsed().as_secs_f64();
        black_box(checksum);
        eprintln!(
            "packed-seq + fxhash:\t{:.2} GB/s",
            LEN as f64 / 1e9 / elapsed
        );
    }

    // nthash crate
    {
        let t = Instant::now();
        let checksum: u64 = nthash::NtHashForwardIterator::new(&ascii, K)
            .unwrap()
            .fold(0u64, |acc, h| acc.wrapping_add(h));
        let elapsed = t.elapsed().as_secs_f64();
        black_box(checksum);
        eprintln!("nthash crate:\t\t{:.2} GB/s", LEN as f64 / 1e9 / elapsed);
    }

    // seq-hash nthash
    {
        let hasher = NtHasher::<false, 1>::new(K);
        let simd_len = (LEN + 7 * (K - 1)).div_ceil(8);
        let t = Instant::now();
        let mut checksum = u32x8::ZERO;
        hasher
            .hash_kmers_simd(packed.as_slice(), 1)
            .advance_with(simd_len, |v| {
                checksum += v;
            });
        let checksum = checksum.as_array_ref().iter().sum::<u32>();
        let elapsed = t.elapsed().as_secs_f64();
        black_box(checksum);
        eprintln!("seq-hash nthash:\t{:.2} GB/s", LEN as f64 / 1e9 / elapsed);
    }

    // seq-hash mulhash
    {
        let hasher = MulHasher::<false, 1>::new(K);
        let simd_len = (LEN + 7 * (K - 1)).div_ceil(8);
        let t = Instant::now();
        let mut checksum = u32x8::ZERO;
        hasher
            .hash_kmers_simd(packed.as_slice(), 1)
            .advance_with(simd_len, |v| {
                checksum += v;
            });
        let checksum = checksum.as_array_ref().iter().sum::<u32>();
        let elapsed = t.elapsed().as_secs_f64();
        black_box(checksum);
        eprintln!("seq-hash mulhash:\t{:.2} GB/s", LEN as f64 / 1e9 / elapsed);
    }
}
