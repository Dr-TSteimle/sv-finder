# sv-finder

```
Usage: sv-finder [OPTIONS] --bam-path <BAM_PATH> --fasta-ref-path <FASTA_REF_PATH> --cytobands-path <CYTOBANDS_PATH>

Options:
  -b, --bam-path <BAM_PATH>                 
  -f, --fasta-ref-path <FASTA_REF_PATH>     
  -o, --output-prefix <OUTPUT_PREFIX>       output file(s) prefix [default: sv-finder]
  -p, --threads <THREADS>                   [default: 1]
  -c, --cytobands-path <CYTOBANDS_PATH>     a cytoband decompressed file, tsv with columns: contig, start, end, cytoband, content.
                                            	- hg19: https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBandIdeo.txt.gz
                                            	- hg38: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBandIdeo.txt.gz
      --distance-threshold <DISTANCE_THRE>  maximum distance in nucleotids between two misaligned reads
                                            for trying to assemble them together
                                             [default: 400]
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
  -h, --help                                Print help (see more with '--help')
  -V, --version                             Print version
  ```
