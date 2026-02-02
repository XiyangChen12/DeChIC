# SCENE

SCENE is a single-cell DNA Deaminase based chromatin immuno-conversion sequencing data analysis pipeline. 1) Pre-processing: Read filtering and alignment, conversion-rate calculation, and removal of SNP-confounded sites; 2) Bulk-level analysis: Peak calling by comparing DeChIC-seq data targeting histone modification or TF to an IgG control; 3) Single-cell analysis: Peak calling for scDeChIC-seq by comparing each cell (or pseudo-bulk signal) to a pseudo-IgG background. SCENE generates a peak-by-cell matrix derived from conversion rates. The matrix can be analyzed with Signac (Seurat) for dimensionality reduction, clustering, and cell-type annotation. SCENE also supports peak clustering based on the proportion of cells with detectable occupancy for each peak.

![Pipeline Diagram](images/SCENE_workflow.pdf)

---

## Installation
```bash
git clone https://github.com/XiyangChen12/DeChIC.git
cd DeChIC
pip install .
```
---

## Usage

Once installed, you can run:
```bash
SCENE -h
usage: SCENE [-h] [-v] {rateMatrix,callpeak,cellMatrix,peakPresence} ...
SCENE is a single-cell DNA Deaminase based chromatin immuno-conversion sequencing data analysis pipeline

positional arguments:
  {rateMatrix,callpeak,cellMatrix,peakPresence}
                        Available subcommands
    rateMatrix          Pre-processing to generate conversion-rate matrix of covered TC sites for each sample (and pseudo-bulk for single-cell mode).
    callpeak            Call peaks for each sample (and pseudo-bulk in single-cell mode).
    cellMatrix          Build a peak-by-cell matrix from the pseudo-bulk peak set.
    peakPresence        Compute Peak Presence Percentage.

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

### Run rateMatrix
```bash
SCENE rateMatrix [-h] -s {bulk,singlecell} --path DIR --SNP BED --genome_fa FASTA
                        [--min_insert_size <int>] [--max_insert_size <int>]
                        --report_repeat_hits {0,1,2} [--max_mismatch <float>]
                        [--seed_size <int>] [--random_seed <int>]
                        --strand_info {0,1,2} [--max_gap <int>] [--convert <str>]

Required Parameters:
--path : directory containing raw FASTQ files
--SNP : BED file of CT/GA SNP sites (C->T and G->A)
--genome_fa : reference genome FASTA (used by BASAL -d and AvgMod)
-s, --sample_type : bulk or singlecell

Optional Parameters (BASAL align):
--min_insert_size : BASAL -m, minimal insert size allowed
--max_insert_size : BASAL -x, maximal insert size allowed
--report_repeat_hits : BASAL -r, how to report repeat hits (0 unique only; 1 random; 2 all)
--max_mismatch : BASAL -v, max mismatch fraction/number
--seed_size : BASAL -s, seed size
--random_seed : BASAL -S, RNG seed for selecting multiple hits
--strand_info : BASAL -n, mapping strand information
--max_gap : BASAL -g, max gap size (<= 1 bp)
--convert : BASAL -M, convert-from:convert-to (e.g. "C:T")
```

### Run callpeak
```bash

usage: SCENE callpeak [-h] -s {bulk,singlecell} --path DIR [-IgG BEDGRAPH]
                      [--met_M <n>] [--met_m <n>] [--met_t <n>] [--p_th <float>]
                      [--d_th <float>] [--merge_d <int>] [--bl_f <float>]
If -s singlecell, this step also calls peaks for pseudo-bulk

Required Parameters:
--path : directory containing bedgraph files
-s, --sample_type : bulk or singlecell

Optional Parameters:
-IgG : IgG bedgraph path. If omitted, SCENE will use pseudo-IgG (mean background) for each sample.
--met_M <n>           -M, --maxdist <n>      maximum distance (default: 50)
--met_m <n>           -m, --mincpgs <n>      minimum CpGs (default: 6)
--met_t <n>           -t, --threads <n>      number of threads (default: 40)
--p_th <float>        Maximum p-value (<) for output peaks (default: 0.05)
--d_th <float>        Minimum mean methylation difference (>=) (default: 0.2)
--merge_d <int>       Minimum distance for merging nearby peaks (default: 200)
--bl_f <float>        Minimum overlap fraction with blacklist; peaks exceeding this fraction will be removed (default: 0.5)
Notes:
  - The options starting with "met_" are passed to the metilene program.
  - filter/merge/blacklist parameters control peak filtering and merging.
```

### Run cellMatrix
```bash
usage: SCENE cellMatrix [-h] --path DIR --peakset BED --output FILE

Required Parameters:
--path : directory containing per-cell bedgraphs
--peakset :pseudo-bulk peak BED file 
--output : output peak-by-cell matrix path
```

### Run peakPresence
```bash
usage: SCENE peakPresence [-h] --path FILE --output FILE

Required Parameters:
--path : peak-by-cell matrix file
--output : output matrix file
```

## Dependencies
External tools required in your environment:

trim_galore
BASAL
samtools
BEDTools
metilene
Python

## Contact
For any questions or suggestions, please contact:
Author: Xiyang Chen
Email: chenxiyang12@outlook.com
GitHub: XiyangChen12
