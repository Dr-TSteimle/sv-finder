
![example workflow](https://github.com/Dr-TSteimle/sv-finder/actions/workflows/rust.yml/badge.svg)

# sv-finder :mag_right: :dna:

Software for genomic structural variant detection and description.

## Quickstart
### :sparkles: Installation

1. Download the latest release:
```
wget https://github.com/Dr-TSteimle/sv-finder/releases/download/1.1.0/sv-finder
```
2. Make the file executable:
```
chmod +x sv-finder
```

You can also install the latest version directly from the GitHub repository using `cargo`:

```bash
cargo install --git https://github.com/Dr-TSteimle/sv-finder
```

### :fire: Usage
```
./sv-finder -h
```

```
Usage: sv-finder [OPTIONS] --bam-path <BAM_PATH> --fasta-ref-path <FASTA_REF_PATH> --cytobands-path <CYTOBANDS_PATH>

Options:
  -b, --bam-path <BAM_PATH>                 
  -f, --fasta-ref-path <FASTA_REF_PATH>     
  -o, --output-prefix <OUTPUT_PREFIX>       output file prefix [default: sv-finder]
  -p, --threads <THREADS>                   [default: 1]
  -c, --cytobands-path <CYTOBANDS_PATH>     a cytoband decompressed file, tsv with columns: contig, start, end, cytoband, content.
                                            	- hg19: https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBandIdeo.txt.gz
                                            	- hg38: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBandIdeo.txt.gz
      --distance-threshold <DISTANCE_THRE>  maximum distance in nucleotids between two misaligned reads
                                            for trying to assemble them together
                                             [default: 350]
      --min-overlapping <MIN_OVERLAPPING>   minimum overlapping length (nt) required for assembling
                                            two reads together
                                             [default: 50]
      --max-consecutive <MAX_CONSECUTIVE>   maximum number of consecutive overlapping mismatch
                                            allowed for assembling reads together
                                             [default: 1]
      --max-mismatches <MAX_DIFFS>          maximum number of overlapping mismatch allowed for
                                            assembling reads together
                                             [default: 3]
      --min-reads <MIN_READS>               minimum reads per cluster
                                             [default: 10]

  -r, --repeat-masker-path <REPEAT_MASKER>
      repeat masker gz file path from UCSC (rmsk.txt.gz)
            - hg19: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
            - hg38: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
  -h, --help                                Print help (see more with '--help')
  -V, --version                             Print version
  ```
sv-finder minimal needed parameters are :
* A sorted BAM file with it's corresponding BAI index file in the same directory.
* The reference FASTA file used for alignement.
* A cytoband TSV file with coordinates from the same reference file (your should decompress and then recompress the file with bgzip).
    * For hg19 you can download it from: https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBandIdeo.txt.gz
    * And for hg38: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBandIdeo.txt.gz

Defaults parameters should be adequate for small reads analysis.

The algorithm is multi-threaded and written in rust :crab: for faster processing.

## :key: Results

The output is a TSV file with the following columns:
* Name of the cluster.
* Genomic coordinates relatives to the reference and to the query (contig:reference_start-reference_stop|query_start-query_stop), the ranges are 1-based and inclusive.
* The [hgvs](https://varnomen.hgvs.org/recommendations/DNA/variant/complex/) unique variant identifier.
* The assembled sequence.

Exemple of a t(10;14) TLX1-TRD observed in an ALL-T sample with it's derivative:
```
chr10:821|chr14:80_0	chr14:22907700-22908008|1-309;chr10:102890895-102890934|320-359	chr14:g.22908008_qter[CGTAGCCCCC;chr10:g.102890895_qter]	GTTAATACTTTACAGTTTTATTACTAGAGGGTTAAAATCCTTTTTCAAGTCTGATAATCAATGATTAACTTTCTTCATTTGTCCTTCACCCATTTGTTTTTTAGGTTGATGGTGTTTTACTTATTGATTTGTGTAATTATAATAATTTTGTGTCTGAGTTTTACAGCATTTAACCACAAAAACAGCATTGGTGAAAGGAGTTTCAGGGGTATTGTGGATGGCAGCGGGTGGTGATGGCAAAGTGCCAAGGAAAGGGAAAAAGGAAGAAGAGGGTTTTTATACTGATGTGTTTCATTGTGCCTTCCTACCGTAGCCCCCGATCTCTGGCTCCGGCATCTGTCTCGGCTTCTGGCGTTCCTGGCCCGCGCGGCGGGCCGCCCTC
chr10:821|chr14:81_1	chr14:22918338-22918105|1-234;chr10:102890891-102890827|240-304	chr14:g.22918105_qterinv[TACCG;chr10:g.102890891_qterinv]	ACCCAAGGAAGAACAGCAGTGAGTGAGAGGTCAGCAGCTGTGGTCATCTCCCTGGTCCAGTCAACTTCCTGCTATCCCTTCCAGGCCCCAAAGCAGGGAGGGAAGCTGCTTGCTGTGTTTGTCTCCTGAGGCATGGGACCCAGGGTGAGGATATCCCAGGGAAATGGCACTTTTGCCCCTGCAGTTTTTGTACAGGTCTCTGTAGGTTTTGTAGCACTGTGCGTATCCCCCAGTACCGTGGGACGGAGACCAAGACTCGGAGTAGTTCATGAAGAGAGAGAAGAGGGGAACAAGGCGAGGCTTA
```

