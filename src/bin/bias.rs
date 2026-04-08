#[path = "../simd.rs"]
mod simd;

use clap::Parser;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use rayon::{ThreadPoolBuilder, current_num_threads};
use serde::Serialize;

use seq_hash::packed_seq::{PackedSeqVec, SeqVec};
use seq_hash::{KmerHasher, NtHasher};

use core::cell::RefCell;
use std::time::Instant;

const DEFAULT_LEN: usize = 500_000_000;
const HIST_SZ: usize = u32::BITS as usize + 1;

thread_local! {
    static CACHE: RefCell<Vec<u32>> = {
        RefCell::new(Vec::with_capacity(DEFAULT_LEN))
    };
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (FASTA, possibly compressed) [default: random seq]
    #[arg(short, long)]
    input: Option<String>,
    /// Output file for histogram (JSON) [default: stdout]
    #[arg(short, long)]
    output: Option<String>,
    /// NtHash R parameter (number of rotations)
    #[arg(short, num_args = 1.., default_value = "1")]
    r: Vec<usize>,
    /// K-mer size
    #[arg(short, num_args = 1.., default_value = "31")]
    k: Vec<usize>,
    /// Hash seed [default: none]
    #[arg(short, long)]
    seed: Option<u32>,
    /// Number of threads [default: all]
    #[arg(short, long)]
    threads: Option<usize>,
}

struct Hist([[usize; HIST_SZ]; HIST_SZ]);

impl Default for Hist {
    fn default() -> Self {
        Self([[0; HIST_SZ]; HIST_SZ])
    }
}

impl std::ops::AddAssign for Hist {
    fn add_assign(&mut self, rhs: Self) {
        for (row, rhs_row) in self.0.iter_mut().zip(rhs.0.iter()) {
            for (a, b) in row.iter_mut().zip(rhs_row.iter()) {
                *a += b;
            }
        }
    }
}

impl Serialize for Hist {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeSeq;
        let mut seq = s.serialize_seq(Some(self.0.len()))?;
        for row in &self.0 {
            seq.serialize_element(row.as_slice())?;
        }
        seq.end()
    }
}

#[derive(Serialize)]
struct HistEntry {
    r: usize,
    k: usize,
    seed: Option<u32>,
    hist: Hist,
}

fn build_hasher(k: usize, r: usize, seed: Option<u32>) -> Box<dyn HistHasher> {
    macro_rules! make {
        ($r:literal) => {
            match seed {
                Some(s) => {
                    Box::new(NtHasher::<false, $r>::new_with_seed(k, s)) as Box<dyn HistHasher>
                }
                None => Box::new(NtHasher::<false, $r>::new(k)),
            }
        };
    }
    match r {
        1 => make!(1),
        2 => make!(2),
        3 => make!(3),
        4 => make!(4),
        5 => make!(5),
        6 => make!(6),
        7 => make!(7),
        8 => make!(8),
        9 => make!(9),
        10 => make!(10),
        11 => make!(11),
        12 => make!(12),
        13 => make!(13),
        14 => make!(14),
        15 => make!(15),
        16 => make!(16),
        _ => panic!("unsupported R={r}, must be in 1..=16"),
    }
}

trait HistHasher {
    fn compute_hist(&self, seq: &PackedSeqVec) -> Hist;
}

impl<H: KmerHasher> HistHasher for H {
    fn compute_hist(&self, seq: &PackedSeqVec) -> Hist {
        let hash_iter = self.hash_kmers_simd(seq.as_slice(), 1);
        let mut hist = Hist::default();
        CACHE.with_borrow_mut(|lz_buf| {
            lz_buf.clear();
            hash_iter.map(simd::lzcnt_u32x8).collect_into(lz_buf);
            lz_buf
                .iter()
                .zip(lz_buf.iter().skip(1))
                .for_each(|(&i, &j)| {
                    hist.0[i as usize][j as usize] += 1;
                });
        });
        hist
    }
}

fn compute_hist(seqs: &[PackedSeqVec], k: usize, seed: Option<u32>, r: usize) -> Hist {
    let mut hists = Vec::with_capacity(seqs.len());
    seqs.par_iter()
        .map(|seq| build_hasher(k, r, seed).compute_hist(seq))
        .collect_into_vec(&mut hists);
    hists.into_iter().fold(Hist::default(), |mut acc, h| {
        acc += h;
        acc
    })
}

fn main() {
    let args = Args::parse();
    let threads = if let Some(t) = args.threads {
        ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .unwrap();
        t
    } else {
        current_num_threads()
    };

    let seqs = if let Some(path) = args.input {
        let mut seqs = vec![PackedSeqVec::default(); threads];
        let mut i = 0;
        let mut reader = parse_fastx_file(path).unwrap();
        while let Some(record) = reader.next() {
            seqs[i].push_ascii(&record.unwrap().seq());
            i += 1;
            if i == threads {
                i = 0;
            }
        }
        seqs
    } else {
        (0..threads)
            .map(|_| PackedSeqVec::random(DEFAULT_LEN))
            .collect()
    };

    let input_size = seqs.iter().map(|seq| seq.len()).sum::<usize>();
    let start = Instant::now();

    let mut res = Vec::with_capacity(args.r.len() * args.k.len());
    for &r in &args.r {
        for &k in &args.k {
            res.push(HistEntry {
                r,
                k,
                seed: args.seed,
                hist: compute_hist(&seqs, k, args.seed, r),
            });
        }
    }

    let time = start.elapsed().as_secs_f64();
    eprintln!(
        "Computed hashes leading zeros in {time:.2}s - {:.2} GB/s",
        (input_size * args.r.len() * args.k.len()) as f64 / 1e9 / time
    );

    println!("{}", serde_json::to_string(&res).unwrap());
}
