// Copyright (C) 2025 [Thomas Steimlé]
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.


use bam::record::tags::TagValue::String as BamString;
use desc_sequence::*;
use noodles_bgzf as bgzf;

use merge_sequences::*;
use rayon::prelude::*;
use rust_lapper::{Interval, Lapper};
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{self, BufRead, BufReader, Write},
    time::Instant,
};

extern crate rand;
use rand::{seq::SliceRandom, rng};

use clap::Parser;
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    bam_path: String,

    #[arg(short, long)]
    fasta_ref_path: String,

    #[arg(
        short,
        long,
        default_value = "sv-finder",
        long_help = "output file(s) prefix"
    )]
    output_prefix: String,

    #[arg(short = 'p', long, default_value_t = 1)]
    threads: usize,

    #[arg(
        short,
        long,
        long_help = "a cytoband decompressed file, tsv with columns: contig, start, end, cytoband, content.\n\t- hg19: https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBandIdeo.txt.gz\n\t- hg38: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBandIdeo.txt.gz"
    )]
    cytobands_path: String,

    #[arg(
        long="distance-threshold",
        default_value_t = 350,
        value_parser = clap::value_parser!(i32).range(1..),
        long_help = "maximum distance in nucleotids between two misaligned reads\nfor trying to assemble them together\n"
    )]
    distance_thre: i32,

    #[arg(
        long,
        long_help="minimum overlapping length (nt) required for assembling\ntwo reads together\n",
        default_value_t = 50,
        value_parser = clap::value_parser!(i16).range(10..),
    )]
    min_overlapping: i16,

    #[arg(
        long, default_value_t = 1,
        long_help="maximum number of consecutive overlapping mismatch\nallowed for assembling reads together\n",
        value_parser = clap::value_parser!(i16).range(1..),
    )]
    max_consecutive: i16,

    #[arg(
        long="max-mismatches", default_value_t = 3,
        long_help="maximum number of overlapping mismatch allowed for\nassembling reads together\n",
        value_parser = clap::value_parser!(i8).range(1..),
    )]
    max_diffs: i8,

    #[arg(
        long="min-reads", default_value_t = 10,
        long_help="minimum reads per cluster\n",
        value_parser = clap::value_parser!(i16).range(1..),
    )]
    min_reads: i16,

    #[arg(
        short = 'r',
        long = "repeat-masker-path",
        long_help = "repeat masker gz file path from UCSC (rmsk.txt.gz)\n\t- hg19: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz\n\t- hg38: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
    )]
    repeat_masker: Option<String>,
}
// target/debug/sv-finder -b /Turbine-pool/LAL-T/data/V4-NextSeq/200604-SureSelect_Capture_V4/2005N140862-BETYA/aln.sorted.bam -f /home/thomas/NGS/ref/hg19.fa -p 30 -c cytoBandIdeo_hg19.txt -o new
//  new_result.txt | grep "chr2:g.pter_43454037[chr14:g.22918105_qterinv]"

fn main() {
    let now = Instant::now();

    // params
    let args = Args::parse();

    let bam_path = &args.bam_path;
    let distance_thre = args.distance_thre;

    let flags: Vec<u16> = vec![81, 161, 97, 145, 65, 129, 113, 177];

    // init the bam reader
    println!("[Bam parser] parsing bam file: {}", bam_path);
    let mut reader_builder = bam::bam_reader::IndexedReaderBuilder::new();
    let parser_threads = if args.threads >= 2 { 2 } else { 1 };
    reader_builder.additional_threads(parser_threads);

    let mut reader = bam::IndexedReader::from_path(bam_path).unwrap();

    // Filter reads with flag 81, 161, 97, 145, 65, 129, 113 or 177
    // or with Supplementary Alignement
    let reader = reader.full_by(move |record| {
        let f = record.flag();
        if flags.contains(&f.0) {
            return true;
        }
        if record.tags().get(b"SA").is_some() {
            return true;
        }
        false
    });

    // Parse the bam header and store correspondance between contigs ids and names
    let header = reader.header();
    let mut id2name: HashMap<i32, String> = HashMap::new();
    for ref_name in header.reference_names().iter() {
        let id = header.reference_id(ref_name).unwrap();
        id2name.insert(id as i32, ref_name.to_string());
    }

    // initialize an HashMap with a entry for each contig, store reads names
    let mut contigs_pos: HashMap<String, Vec<(String, i32)>> = HashMap::new();
    id2name.iter().for_each(|(_, contig)| {
        contigs_pos.insert(contig.to_string(), Vec::new());
    });

    let mut abn_reads: HashSet<String> = HashSet::new();
    let mut sequences: HashMap<String, String> = HashMap::new();
    let mut n_abn_positions = 0;

    reader.for_each(|result| {
        if let Ok(record) = result {
            let qname = String::from_utf8_lossy(record.name()).to_string();
            // let read_sens = if record.flag().first_in_pair() { "_R1" } else { "_R2" };
            let key_base = qname.clone();
            let mut iter_key = 0;

            // add to sequences
            loop {
                if sequences.contains_key(&format!("{}_{}", key_base, iter_key)) {
                    iter_key += 1;
                } else {
                    break;
                }
            }
            let k = format!("{}_{}", key_base, iter_key);

            if !sequences.contains_key(&k) {
                let mut sequence = record.sequence().to_vec();
                if record.flag().is_reverse_strand() {
                    sequence = record
                        .sequence()
                        .rev_compl(0..sequence.len())
                        .collect::<Vec<u8>>();
                }
                sequences.insert(k.clone(), String::from_utf8_lossy(&sequence).to_string());
            }

            let pos = record.start();
            let contig = id2name.get(&record.ref_id()).unwrap();

            // add to contigs_pos
            if abn_reads.insert(format!("{} {}:{}", qname, contig, pos)) {
                let contig_pos = contigs_pos.get_mut(contig).unwrap();
                contig_pos.push((qname.to_string(), pos));
            };
            n_abn_positions += 1;

            // SA
            if let Some(BamString(sa, _)) = record.tags().get(b"SA") {
                let (contig, pos) = process_sa(sa);
                if abn_reads.insert(format!("{} {}:{}", qname, contig, pos)) {
                    contigs_pos
                        .entry(contig)
                        .or_default()
                        .push((qname.to_owned(), pos));
                    n_abn_positions += 1;
                };
            }

            if (n_abn_positions + 1) % 1_000 == 0 {
                println!(
                    "[Bam parser] {} abnormal reads parsed in {:#?}",
                    n_abn_positions + 1,
                    now.elapsed()
                );
            }
        }
    });
    println!(
        "[Bam parser] {} abnormal reads parsed in {:#?}",
        n_abn_positions,
        now.elapsed()
    );

    // Clustering the reads based on their positions,
    // reads with the same set of positions +/- distance_threshold
    let mut grps = contigs_pos.clone();
    for (contig, pos) in contigs_pos.iter_mut() {
        if !pos.is_empty() {
            pos.sort_by(|a, b| a.1.cmp(&b.1));

            let mut last_group = 0;
            let mut last_pos = pos[0].1;

            let grps_loc = grps.get_mut(contig).unwrap();

            *grps_loc = pos
                .iter()
                .map(|(qname, pos)| {
                    if pos - last_pos > distance_thre {
                        last_group += 1;
                    }
                    last_pos = *pos;
                    (qname.to_string(), last_group)
                })
                .collect::<Vec<(String, i32)>>();
        }
    }

    let mut by_qnames: HashMap<String, Vec<String>> = HashMap::new();
    for (contig, grps) in grps.iter() {
        for (qname, grp) in grps.iter() {
            let v = format!("{}:{}", contig, grp);
            if let Some(qnames_grps) = by_qnames.get_mut(qname) {
                qnames_grps.push(v);
                qnames_grps.dedup();
            } else {
                by_qnames.insert(qname.to_string(), vec![v]);
            }
        }
    }

    let mut by_clust: HashMap<String, Vec<String>> = HashMap::new();
    for (qname, cluster) in by_qnames.iter_mut() {
        cluster.sort();
        let k = cluster.join("|");
        if let Some(qnames) = by_clust.get_mut(&k) {
            qnames.push(qname.to_string());
        } else {
            by_clust.insert(k, vec![qname.to_string()]);
        }
    }

    let min_reads = args.min_reads as usize;
    let mut clusters: Vec<(String, Vec<String>)> = Vec::from_iter(by_clust)
        .into_iter()
        .filter(|x| x.1.len() >= min_reads)
        .collect();

    clusters.sort_by(|a, b| b.1.len().cmp(&a.1.len()));

    println!(
        "[Clustering] {} clusters formed in {:#?}",
        clusters.len(),
        now.elapsed()
    );

    // Merging reads for each cluster in parallel
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    let res: Vec<(String, Vec<DNAString>, (usize, usize))> = clusters // let res: Vec<(String, Vec<DNAString>)> = clusters
        .par_iter()
        .enumerate()
        .map(|(i, (cluster_name, reads))| {
            let mut seqs = Vec::new();

            println!("[Merging] cluster {} {}", i, cluster_name);
            for read in reads.iter() {
                let mut iter_key = 0;

                loop {
                    let k = format!("{}_{}", read, iter_key);

                    if sequences.contains_key(&k) {
                        let s = sequences.get(&k).unwrap().to_string();
                        if s.chars().all(|x| ['A', 'T', 'C', 'G'].contains(&x)) {
                            seqs.push(DNAString::new(s.as_bytes().to_vec()));
                        }
                    } else {
                        break;
                    }
                    iter_key += 1;
                }
            }

            let kmer_len = args.min_overlapping as usize;

            seqs.retain(|e| e.0.len() >= kmer_len);
            let total_seqs = seqs.len();
            println!(
                "[Merging] cluster {} {}: {} sequences to merge",
                i, cluster_name, total_seqs
            );

            if seqs.len() > 1 {
                let mut last_len = 0;
                let mut iter = 0;
                while last_len != seqs.len() {
                    seqs = rec_merger(seqs, kmer_len, args.max_diffs, args.max_consecutive as i8);
                    last_len = seqs.len();
                    iter += 1;
                    if iter > 10 {
                        break;
                    }
                }
            }

            println!("[Merging] cluster {} {}: Finished", i, cluster_name);

            let stats = (total_seqs, seqs.len());

            (cluster_name.to_string(), seqs, stats)
        })
        .collect();

    let mut res: Vec<(String, DNAString)> = res
        .into_iter()
        .flat_map(|(name, merged, _)| {
            merged
                .into_iter()
                .enumerate()
                .map(|(i, e)| (format!("{}_{}", name, i), e))
                .collect::<Vec<(String, DNAString)>>()
        })
        .collect();

    res.sort_by(|a, b| b.1 .0.len().cmp(&a.1 .0.len()));
    println!("{} clusters merged in {:#?}", res.len(), now.elapsed());

    // Use bwa mem aligner to describe new contigs and forge an unique hgvs
    let mut result_file = File::create(format!("{}_result.txt", args.output_prefix)).unwrap();
    let centro_pos = parse_cytoband(&args.cytobands_path);
    let aligner = BwaAlign::init(&args.fasta_ref_path);
    let lappers = if let Some(path) = &args.repeat_masker {
        load_repeat(path).unwrap()
    } else {
        HashMap::new()
    };
    res.iter().for_each(|(name, sequence)| {
        let mut res = Result::new(name.to_string(), sequence.as_string());
        // res.align(&args.fasta_ref_path, &centro_pos);
        let ranges = aligner.get_ref_positions(&sequence.as_string());
        res.align_bwa(ranges, &centro_pos);

        if args.repeat_masker.is_some() {
            res.check_repeat(&lappers);
        }
        if res.hgvs.is_some() {
            writeln!(result_file, "{}", res.print_tsv_line()).unwrap();
        }
    });

    println!("Process done in {:#?}", now.elapsed());
}

fn process_sa(sa: &[u8]) -> (String, i32) {
    let mut contig = String::default();
    let mut pos: i32 = i32::default();

    for sa in sa.split(|c| b';' == *c).filter(|&x| !x.is_empty()) {
        for (i, sa_part) in sa
            .split(|c| b',' == *c)
            .filter(|&x| !x.is_empty())
            .enumerate()
        {
            match i {
                0 => contig = String::from_utf8_lossy(sa_part).to_string(),
                1 => {
                    pos = String::from_utf8_lossy(sa_part)
                        .to_string()
                        .parse::<i32>()
                        .unwrap()
                }
                _ => break,
            };
        }
    }
    (contig, pos)
}

pub fn merger(
    sequences: &[DNAString],
    to_merge: &DNAString,
    kmer_len: usize,
    max_mismatches: i8,
    max_consecutive_mismatches: i8,
) -> Option<(DNAString, Vec<DNAString>)> {
    let mut to_merge = to_merge.clone();
    let mut merged = Vec::new();

    // try to merge to others seeds
    for (i, seed) in sequences.iter().enumerate() {
        if seed.0.len() >= kmer_len {
            let r = to_merge.merge(seed, kmer_len, max_mismatches, max_consecutive_mismatches);
            match r {
                Ok(_) => {
                    merged.push(i);
                }
                Err(_err) => (),
            }
        }
    }

    if !merged.is_empty() {
        Some((
            to_merge,
            sequences
                .iter()
                .enumerate()
                .filter(|(i, _)| !merged.contains(i))
                .map(|(_, e)| e.clone())
                .collect::<Vec<DNAString>>(),
        ))
    } else {
        None
    }
}

pub fn rec_merger(
    mut seqs: Vec<DNAString>,
    kmer_len: usize,
    max_mismatches: i8,
    max_consecutive_mismatches: i8,
) -> Vec<DNAString> {
    let mut n_retry = 0;

    let mut merged = Vec::new();
    let mut last_len = seqs.len();
    loop {
        if last_len <= 1 {
            break;
        }

        let tmp = seqs.get(1..).unwrap().to_vec();
        let re = merger(
            &tmp,
            seqs.first().unwrap(),
            kmer_len,
            max_mismatches,
            max_consecutive_mismatches,
        );
        match re {
            Some(res) => {
                merged.push(res.0.clone());
                // res.1.append(&mut vec![res.0]);
                seqs = res.1;
            }
            None => match n_retry {
                0 => seqs.reverse(),
                _ => seqs.shuffle(&mut rng()),
            },
        }
        n_retry += 1;

        if last_len == seqs.len() {
            break;
        } else {
            last_len = seqs.len();
        }
    }
    merged
}

type FullRange = ((String, i32, i32), (i32, i32));
#[derive(Debug, Clone)]
pub struct Result {
    pub name: String,
    pub sequence: String,
    pub ranges: Option<Vec<FullRange>>,
    pub hgvs: Option<Vec<String>>,
    pub on_repeat: Option<bool>,
}

impl Result {
    pub fn new(name: String, sequence: String) -> Result {
        Result {
            name,
            sequence,
            ranges: None,
            hgvs: None,
            on_repeat: None,
        }
    }

    pub fn align_bwa(
        &mut self,
        mut ranges: Vec<FullRange>,
        centro_pos: &[(String, i32)],
    ) {
        ranges.sort_by(|a, b| a.1 .0.cmp(&b.1 .0));

        if !ranges.is_empty() {
            self.ranges = Some(ranges);
            self.hgvs(centro_pos);
        }
    }

    pub fn hgvs(&mut self, centro_pos: &[(String, i32)]) {
        if let Some(ranges) = self.ranges.clone() {
            if ranges.len() >= 2 {
                let mut lag_ranges: Option<&FullRange> = None;
                let mut hgvs: Vec<String> = Vec::new();
                for range in ranges.iter() {
                    if let Some(a_range) = lag_ranges {
                        let (a_contig, a_first, a_pos) = a_range.0.to_owned();
                        let (b_contig, b_first, b_pos) = range.0.to_owned();

                        let a_rc = a_first > a_pos;
                        let b_rc = b_first > b_pos;

                        let ins = substring_inc(self.sequence.clone(), a_range.1 .1, range.1 .0);

                        if let Some(res) = hgvs_delins(
                            a_contig, a_pos, a_rc, b_contig, b_first, b_rc, ins, centro_pos,
                        ) {
                            hgvs.push(res);
                        }

                        lag_ranges = Some(range);
                    } else {
                        lag_ranges = Some(range);
                    }
                }
                if !hgvs.is_empty() {
                    self.hgvs = Some(hgvs);
                }
            }
        }
    }

    pub fn print_tsv_line(&self) -> String {
        let mut str_ranges = Vec::new();
        if let Some(ranges) = self.ranges.clone() {
            for range in ranges {
                str_ranges.push(format!(
                    "{}:{}-{}|{}-{}",
                    range.0 .0, range.0 .1, range.0 .2, range.1 .0, range.1 .1
                ));
            }
        }
        let str_ranges = if !str_ranges.is_empty() {
            str_ranges.join(";")
        } else {
            "NA".to_string()
        };
        let str_hgvs = if let Some(hgvs) = self.hgvs.clone() {
            hgvs.join(";")
        } else {
            "NA".to_string()
        };

        let rep = if let Some(on_rep) = self.on_repeat {
            if on_rep {
                "1"
            } else {
                "0"
            }
        } else {
            "NA"
        };

        [
            self.name.clone(),
            str_ranges,
            str_hgvs,
            self.sequence.clone(),
            rep.to_string(),
        ]
        .join("\t")
    }

    pub fn check_repeat(&mut self, lappers: &HashMap<String, Lapper<u32, String>>) {
        if let Some(ranges) = &self.ranges {
            if ranges.iter().any(|((chr, start, end), _)| {
                if let Some(lapper) = lappers.get(chr) {
                    lapper.find(*start as u32, *end as u32).count() > 0
                } else {
                    false
                }
            }) {
                self.on_repeat = Some(true)
            } else {
                self.on_repeat = Some(false)
            }
        }
    }
}

pub fn parse_cytoband(path: &str) -> Vec<(String, i32)> {
    let lines = io::BufReader::new(File::open(path).unwrap()).lines();
    let mut results: Vec<(String, i32)> = Vec::new();
    for line in lines {
        let line = line.unwrap();
        if line.contains("acen") {
            let mut col = line.split('\t');
            let contig = col.next().unwrap().to_string();
            if results
                .iter()
                .filter(|e| e.0 == contig)
                .collect::<Vec<_>>()
                .is_empty()
            {
                results.push((contig, col.nth(1).unwrap().parse().unwrap()));
            }
        }
    }
    results
}

#[allow(clippy::too_many_arguments)]
pub fn hgvs_delins(
    a_contig: String,
    a_pos: i32,
    a_rc: bool,
    b_contig: String,
    b_pos: i32,
    b_rc: bool,
    ins: Option<String>,
    centro_pos: &[(String, i32)],
) -> Option<String> {
    let a_centro_pos: Vec<i32> = centro_pos
        .iter()
        .filter(|(c, _)| *c == a_contig)
        .map(|(_, p)| *p)
        .collect();
    let b_centro_pos: Vec<i32> = centro_pos
        .iter()
        .filter(|(c, _)| *c == b_contig)
        .map(|(_, p)| *p)
        .collect();

    if a_centro_pos.len() == 1 && b_centro_pos.len() == 1 {
        let a_arm = if a_pos < *a_centro_pos.first().unwrap() {
            "p"
        } else {
            "q"
        };
        let b_arm = if b_pos < *b_centro_pos.first().unwrap() {
            "p"
        } else {
            "q"
        };

        let a_hgvs = hgvs_part(&a_contig, a_pos, a_arm, a_rc);
        let b_hgvs = hgvs_part(&b_contig, b_pos, b_arm, b_rc);

        let ins = match ins {
            Some(i) => format!("{i};"),
            None => "".to_string(),
        };

        Some(format!("{}[{}{}]", a_hgvs, ins, b_hgvs))
    } else {
        None
    }
}

pub fn hgvs_part(contig: &str, pos: i32, arm: &str, rc: bool) -> String {
    match (arm, rc) {
        ("p", false) => format!("{}:g.pter_{}", contig, pos),
        // ("p", true)  => format!("{}:g.{}_pter", contig, pos), // not clear in the hgvs standard
        ("p", true) => format!("{}:g.pter_{}_inv", contig, pos),
        ("q", false) => format!("{}:g.{}_qter", contig, pos),
        ("q", true) => format!("{}:g.{}_qterinv", contig, pos),
        // ("q", true)  => format!("{}:g.qter_{}", contig, pos),
        (&_, _) => "".to_string(),
    }
}

fn substring_inc(seq: String, start: i32, stop: i32) -> Option<String> {
    let len = stop - start - 1;
    if len > 0 && (start - 1) >= 0 {
        Some(
            seq.chars()
                .skip((start - 1) as usize)
                .take(len as usize)
                .collect(),
        )
    } else {
        None
    }
}

fn load_repeat(path: &str) -> anyhow::Result<HashMap<String, Lapper<u32, String>>> {
    // Open the gzipped file
    let file = File::open(path)?;
    let bgzf = bgzf::Reader::new(file);
    let reader = BufReader::new(bgzf);

    // Create a Lapper to store intervals
    let mut hm_lapper: HashMap<String, Lapper<u32, String>> = HashMap::new();

    // Read and parse the file, adding intervals to the Lapper
    let mut curr_lapper: Lapper<u32, String> = Lapper::new(vec![]);
    let mut curr_chrom = String::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 8 {
            let chr: String = parts[5].to_string();
            if curr_chrom.is_empty() {
                println!("[Repeat lapper] {chr}");
                curr_chrom = chr.clone();
            }

            if chr != curr_chrom {
                println!("[Repeat lapper] {chr}");
                hm_lapper.insert(curr_chrom, curr_lapper);
                curr_chrom = chr.clone();
                curr_lapper = Lapper::new(vec![]);
            }
            let start: u32 = parts[6].parse()?;
            let stop: u32 = parts[7].parse()?;
            let name = parts[0].to_string();
            curr_lapper.insert(Interval {
                start,
                stop,
                val: name,
            });
        }
    }

    hm_lapper.insert(curr_chrom, curr_lapper);
    println!("[Repeat lapper] n chr: {}", hm_lapper.len());

    Ok(hm_lapper)
}
