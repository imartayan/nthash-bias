use seq_hash::packed_seq::{AsciiSeqVec, PackedSeqVec, Seq, SeqVec, u32x8};
use seq_hash::{KmerHasher, MulHasher, NtHasher};

use core::hash::BuildHasher;
use core::hint::black_box;
use rustc_hash::FxBuildHasher;
use std::time::Instant;

const LEN: usize = 100_000_000;
const K: usize = 31;
const MASK: u64 = (1 << (2 * K)) - 1;

fn main() {
    let ascii: Vec<u8> = AsciiSeqVec::random(LEN).seq;
    let packed = PackedSeqVec::from_ascii(&ascii);

    // nthash-rs crate
    {
        let t = Instant::now();
        let mut checksum: u64 = 0;
        let mut hasher = nthash_rs::NtHash::new(&ascii, K as u16, 1, 0).unwrap();
        while hasher.roll() {
            checksum = checksum.wrapping_add(hasher.hashes()[0]);
        }
        let elapsed = t.elapsed().as_secs_f64();
        black_box(checksum);
        eprintln!("nthash-rs crate:\t{:.2} GB/s", LEN as f64 / 1e9 / elapsed);
    }

    // nthash crate
    {
        let t = Instant::now();
        let checksum: u64 = nthash::NtHashIterator::new(&ascii, K)
            .unwrap()
            .fold(0u64, |a, h| a.wrapping_add(h));
        let elapsed = t.elapsed().as_secs_f64();
        black_box(checksum);
        eprintln!("nthash crate:\t\t{:.2} GB/s", LEN as f64 / 1e9 / elapsed);
    }

    // ascii iter + fxhash
    {
        let t = Instant::now();
        let checksum = ascii.windows(K).fold(0u64, |acc, kmer| {
            acc.wrapping_add(FxBuildHasher.hash_one(kmer))
        });
        let elapsed = t.elapsed().as_secs_f64();
        black_box(checksum);
        eprintln!(
            "ascii seq + fxhash:\t{:.2} GB/s",
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
                acc.wrapping_add(FxBuildHasher.hash_one(kmer))
            });
        let elapsed = t.elapsed().as_secs_f64();
        black_box(checksum);
        eprintln!(
            "packed-seq + fxhash:\t{:.2} GB/s",
            LEN as f64 / 1e9 / elapsed
        );
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
