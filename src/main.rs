// Copyright (c) <2023> <Thomas Steimlé>

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

use std::{time::Instant, collections::{HashMap, HashSet}, fs::File, io::{self, BufRead, Write}};
use bam::record::tags::TagValue::String as BamString;
use rayon::prelude::*;
use merge_sequences::*;
use desc_sequence::*;

extern crate rand;
use rand::{thread_rng, seq::SliceRandom};

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
        long_help = "output file prefix"
    )]
    output_prefix: String,
   
    #[arg(short='p', long, default_value_t = 1)]
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
        value_parser = clap::value_parser!(i8).range(1..),
    )]
    min_reads: i8,
}
// target/debug/sv-finder -b /Turbine-pool/LAL-T/data/V4-NextSeq/200604-SureSelect_Capture_V4/2005N140862-BETYA/aln.sorted.bam -f /home/thomas/NGS/ref/hg19.fa -p 30 -c cytoBandIdeo_hg19.txt
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
        if flags.contains(&f.0) { return true }
        if let Some(_) = record.tags().get(b"SA") { return true }
        false
    });

    // Parse the bam header and store correspondance between contigs ids and names
    let header = reader.header();
    let mut id2name: HashMap<i32, String> = HashMap::new();
    for ref_name in header.reference_names().iter() {
        let id = header.reference_id(ref_name).unwrap();
        id2name.insert(id as i32, ref_name.to_string());
    }

    // initialize an HashMap with a entry for each contig that'll store reads names
    let mut contigs_pos: HashMap<String, Vec<(String, i32)>> = HashMap::new();
    id2name.iter().for_each(|(_, contig)| {
        contigs_pos.insert(contig.to_string(), Vec::new());
    });
    
    let mut abn_reads: HashSet<String> = HashSet::new();
    let mut sequences: HashMap<String, String> = HashMap::new();
    let mut n_abn_positions = 0;
    
    for result in reader {
        if let Ok(record) = result {
            let qname    = String::from_utf8_lossy(record.name()).to_string();
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
                    sequence = record.sequence().rev_compl(0..sequence.len()).collect::<Vec<u8>>();
                }
                sequences.insert(k.clone(), String::from_utf8_lossy(&sequence).to_string());
            }
            
            let pos        = record.start();
            let contig = id2name.get(&record.ref_id()).unwrap();
  
            // add to contigs_pos
            if abn_reads.insert(
                format!("{} {}:{}", qname.to_string(), contig.to_string(), pos)
            ) {
                let contig_pos = contigs_pos.get_mut(contig).unwrap();
                contig_pos.push((qname.to_string(), pos));
            };
            n_abn_positions += 1;
            
            // SA
            if let Some(sa) = record.tags().get(b"SA") {
                match sa {
                    BamString(sa, _) => {
                        let (contig, pos) = process_sa(sa);
                        if abn_reads.insert(
                            format!("{} {}:{}", qname.to_string(), contig.to_string(), pos)
                        ) {
                            let contig_pos = contigs_pos.get_mut(&contig).unwrap();
                            contig_pos.push((qname, pos));
                            n_abn_positions += 1;
                        };
                    },
                    _ => {},
                }
            }

            if (n_abn_positions + 1) % 1_000 == 0 {
                println!("[Bam parser] {} abnormal reads parsed in {:#?}", n_abn_positions + 1,  now.elapsed());
            }
        }
    }
    println!("[Bam parser] {} abnormal reads parsed in {:#?}", n_abn_positions,  now.elapsed());

    // Clustering the reads based on their positions,
    // reads with the same set of positions +/- distance_threshold
    let mut grps = contigs_pos.clone();
    for (contig, pos) in contigs_pos.iter_mut() {
        if pos.len() > 0 {
            pos.sort_by(|a,b| a.1.cmp(&b.1));

            let mut last_group = 0;
            let mut last_pos = pos[0].1;
    
            let grps_loc = grps.get_mut(contig).unwrap();
            
            *grps_loc = pos.iter()
            .map(|(qname, pos)| {
                if pos - last_pos > distance_thre {
                    last_group += 1;
                }
                last_pos = *pos;
                (qname.to_string(),last_group)
            }).collect::<Vec<(String, i32)>>();
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
    let mut clusters: Vec<(String, Vec<String>)> = Vec::from_iter(by_clust.into_iter())
        .into_iter()
        .filter(|x| x.1.len() >= min_reads)
        .collect();
    
    clusters
        .sort_by(|a,b| b.1.len().cmp(&a.1.len()));

    println!("[Clustering] {} clusters formed in {:#?}", clusters.len(), now.elapsed());
    
    // Merging reads for each cluster in parallel
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global().unwrap();

    let res: Vec<(String, Vec<DNAString>, (usize, usize))> = clusters// let res: Vec<(String, Vec<DNAString>)> = clusters
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
                        if s.chars().all(|x| vec!['A', 'T', 'C', 'G'].contains(&x)) {
                            seqs.push(DNAString::new(s.as_bytes().to_vec()));
                        }
                    } else {
                        break;
                    }
                    iter_key += 1;
                }
            }
            
            let kmer_len = args.min_overlapping as usize;

            seqs = seqs.into_iter().filter(|e| e.0.len() >= kmer_len).collect();
            let total_seqs = seqs.len();
            println!("[Merging] cluster {} {}: {} sequences to merge", i, cluster_name, total_seqs);

            if seqs.len() > 1 {
                let mut last_len = 0;
                let mut iter = 0;
                while last_len != seqs.len() {
                    seqs = rec_merger(seqs, kmer_len, args.max_diffs, args.max_consecutive as i8);
                    last_len = seqs.len();
                    iter += 1;
                    if iter > 10 { break; }
                }
            }
            
            println!("[Merging] cluster {} {}: Finished", i, cluster_name);

            let stats = (total_seqs, seqs.len());
            
            (cluster_name.to_string(), seqs, stats)
        }).collect();
    
    let mut res: Vec<(String, DNAString)> = res
        .into_iter()
        .flat_map(|(name, merged, _)| {
            merged
                .into_iter()
                .enumerate()
                .map(|(i,e)| {(format!("{}_{}", name, i), e)} )
                .collect::<Vec<(String, DNAString)>>()
        })
        .collect();
    
    res.sort_by(|a, b| {
        b.1.0.len().cmp(&a.1.0.len())
    });
    println!("{} clusters merged in {:#?}", res.len(), now.elapsed());

    // Use bwa mem aligner to describe new contigs and forge an unique hgvs
    let mut result_file = File::create(format!("{}_result.txt", args.output_prefix)).unwrap();
    let centro_pos = parse_cytoband(&args.cytobands_path);
    let aligner = BwaAlign::init(&args.fasta_ref_path);
    res
        .iter()
        .for_each(|(name, sequence)| {
            let mut res = Result::new(name.to_string(), sequence.as_string());
            // res.align(&args.fasta_ref_path, &centro_pos);
            let ranges = aligner.get_ref_positions(&sequence.as_string());
            res.align_bwa(ranges, &centro_pos);
            
            if res.hgvs.is_some() {
                writeln!(result_file, "{}", res.print_tsv_line()).unwrap();
            } else {
                let ranges = if let Some(ranges) = res.ranges {
                    let ranges: Vec<String> = ranges.iter().map(|((c, rs, re), (ss, se))| format!("{}:{}-{}|{}-{}", c, rs, re, ss, se)).collect();
                    ranges.join(";")
                } else {
                    "".to_string()
                };
                writeln!(result_file, "{}\t{}\t\t{}", name, ranges, sequence.as_string()).unwrap();
            }
        });

    println!("Process done in {:#?}", now.elapsed());
}

fn process_sa (sa: &[u8]) -> (String, i32) {
    let mut contig = String::default();
    let mut pos: i32 = i32::default();

    for sa in sa.split(|c| b';' == *c ).filter(|&x| !x.is_empty()) {
        for (i, sa_part) in sa.split(|c| b',' == *c ).filter(|&x| !x.is_empty()).enumerate() {
            match i {
                0 => contig = String::from_utf8_lossy(sa_part).to_string(),
                1 => pos = String::from_utf8_lossy(sa_part).to_string().parse::<i32>().unwrap(),
                _ => break
            };
        }
    }
    (contig, pos)
}

pub fn merger(sequences: &Vec<DNAString>, to_merge: &DNAString, kmer_len: usize, max_mismatches: i8, max_consecutive_mismatches: i8) -> Option<(DNAString, Vec<DNAString>)> {
    let mut to_merge = to_merge.clone();
    let mut merged = Vec::new();

    // try to merge to others seeds
    for (i, seed) in sequences.iter().enumerate() {
        if seed.0.len() >= kmer_len {
            let r = to_merge.merge(seed, kmer_len, max_mismatches, max_consecutive_mismatches);
            match r {
                Ok(_) => {
                    merged.push(i);
                },
                Err(_err) => (),
            }
        }
    }
    
    if merged.len() > 0 {
        Some((
            to_merge,
            sequences.into_iter().enumerate().filter(|(i, _)| !merged.contains(i)).map(|(_, e)| e.clone()).collect::<Vec<DNAString>>()
        ))
    } else {
        None
    }
}

pub fn rec_merger(mut seqs: Vec<DNAString>, kmer_len: usize, max_mismatches: i8, max_consecutive_mismatches: i8) -> Vec<DNAString> {
    let mut n_retry = 0;

    let mut merged = Vec::new();
    let mut last_len = seqs.len();
    loop {
        if last_len <= 1 { break; }

        let tmp = seqs.get(1..).unwrap().to_vec();
        let re = merger(&tmp, seqs.get(0).unwrap(), kmer_len, max_mismatches, max_consecutive_mismatches);
        match re {
            Some(res) => {
                merged.push(res.0.clone());
                // res.1.append(&mut vec![res.0]);
                seqs = res.1;
            },
            None => {
                match n_retry {
                    0 => seqs.reverse(),
                    _ => seqs.shuffle(&mut thread_rng())
                }
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

#[derive(Debug, Clone)]
pub struct Result {
    pub name    : String,
    pub sequence: String,
    pub ranges  : Option<Vec<((String, i32, i32), (i32, i32))>>,
    pub hgvs    : Option<Vec<String>>
}

impl Result {
    pub fn new(name: String, sequence: String) -> Result {
        Result { name, sequence, ranges: None, hgvs: None }
    }

    pub fn align_bwa(&mut self, mut ranges: Vec<((String, i32, i32), (i32, i32))>, centro_pos: &Vec<(String, i32)>) {
        ranges
        .sort_by(|a, b| a.1.0.cmp(&b.1.0));

        if ranges.len() > 0 {
            self.ranges = Some(ranges);
            self.hgvs(centro_pos);
        }
    }

    pub fn hgvs(&mut self, centro_pos: &Vec<(String, i32)>) {
        if let Some(ranges) = self.ranges.clone() {
            if ranges.len() >= 2 {
                let mut lag_ranges: Option<&((String, i32, i32), (i32, i32))> = None;
                let mut hgvs: Vec<String> = Vec::new();
                for range in ranges.iter() {
                    if let Some(a_range) = lag_ranges {
                        let (a_contig, a_first, a_pos) = a_range.0.to_owned();
                        let (b_contig, b_first, b_pos) = range.0.to_owned();
                        
                        let a_rc = a_first > a_pos;
                        let b_rc = b_first > b_pos;
                        
                        let ins = substring_inc(self.sequence.clone(), a_range.1.1, range.1.0);

                        if let Some(res) =hgvs_delins(a_contig, a_pos, a_rc, b_contig, b_first, b_rc, ins, centro_pos) {
                            hgvs.push(res);
                        }
            
                        lag_ranges = Some(range);
                    } else {
                        lag_ranges = Some(range);
                    }
                }
                if hgvs.len() > 0 {
                    self.hgvs = Some(hgvs);
                }
            }
        }
    }
    
    pub fn print_tsv_line(&self) -> String {
        let mut str_ranges = Vec::new();
        if let Some(ranges) = self.ranges.clone() {
            for range in ranges {
                str_ranges.push(format!("{}:{}-{}|{}-{}", range.0.0, range.0.1, range.0.2, range.1.0, range.1.1));
            }
        }
        let str_ranges = if str_ranges.len() > 0 { str_ranges.join(";") } else { "NA".to_string() };
        let str_hgvs = if let Some(hgvs) = self.hgvs.clone() { hgvs.join(";") } else { "NA".to_string() };

        vec![
            self.name.clone(),
            str_ranges,
            str_hgvs,
            self.sequence.clone()
        ].join("\t")
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
            if results.iter().filter(|e| e.0.to_string() == contig).collect::<Vec<_>>().len() == 0 {
                results.push((contig, col.nth(1).unwrap().parse().unwrap()));
            }
        }
    }
    results
}

pub fn hgvs_delins(
    a_contig: String, a_pos: i32, a_rc: bool,
    b_contig: String, b_pos: i32, b_rc: bool,
    ins     : Option<String>,
    centro_pos: &Vec<(String, i32)>
) -> Option<String> {
    let a_centro_pos: Vec<i32> = centro_pos.iter().filter(|(c, _)| c.to_string() == a_contig).map(|(_, p)| *p).collect();
    let b_centro_pos: Vec<i32> = centro_pos.iter().filter(|(c, _)| c.to_string() == b_contig).map(|(_, p)| *p).collect();

    if a_centro_pos.len() == 1 && b_centro_pos.len() == 1 {
        let a_arm = if a_pos < a_centro_pos.get(0).unwrap().to_owned() { "p" } else { "q" };
        let b_arm = if b_pos < b_centro_pos.get(0).unwrap().to_owned() { "p" } else { "q" };

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
        ("p", true)  => format!("{}:g.pter_{}_inv", contig, pos),
        ("q", false) => format!("{}:g.{}_qter", contig, pos),
        ("q", true)  => format!("{}:g.{}_qterinv", contig, pos),
        // ("q", true)  => format!("{}:g.qter_{}", contig, pos),
        (&_, _) => "".to_string(),
    }
}

fn substring_inc(seq: String, start: i32, stop: i32) -> Option<String> {
    let len = stop - start - 1;
    if len > 0 && (start - 1) >= 0 {
        Some(seq.chars().skip((start - 1) as usize).take(len as usize).collect())
    } else {
        None
    }
}
