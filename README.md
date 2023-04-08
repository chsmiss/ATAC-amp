# ATAC-amp
Searching for  co-amplified regions on the genome from ATAC-seq data
python find_spreads_process_v2.py -h
usage: find_spreads_process_v2.py [-h] [--bam BAM] [--name NAME]
                                  [--isize_value ISIZE]
                                  [--interval_size INTERVAL] [--mapq MAQP]
                                  [--mode {0,1,2}] [--discbk DISCBK]
                                  [--type LIB] [--gtf GTF]

use atac data find ecdna

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM             input the bam file
  --name NAME, -n NAME  prefix of output files
  --isize_value ISIZE, -i ISIZE
                        judge a pair of reads whether is discordant
  --interval_size INTERVAL, -s INTERVAL
                        size of interval when compute breakpoint nearby
                        coverage
  --mapq MAQP, -q MAQP  reads maqp threshold
  --mode {0,1,2}, -m {0,1,2}
                        choose the analysis mode,0/1/2
  --discbk DISCBK, -d DISCBK
                        if you choose mode 1,you need to input discbk file
  --type LIB            choose library:sc/bulk
  --gtf GTF             gtf file
