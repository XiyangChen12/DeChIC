# modules/cellMatrix.snakefile
from pathlib import Path
import os, glob

BASE = Path(workflow.basedir)
PKG  = BASE.parent
configfile: str(PKG / "config.yaml")

# ---- required inputs from CLI ----
ratio_dir = config.get("path", None)       # --path
peakset   = config.get("peakset", None)    # --peakset
out_matrix= config.get("output", None)     # --output

if ratio_dir is None or str(ratio_dir).strip() == "":
    raise ValueError("Missing required --path (bedgraph directory).")
if peakset is None or str(peakset).strip() == "":
    raise ValueError("Missing required --peakset (merged peak bed).")
if out_matrix is None or str(out_matrix).strip() == "":
    raise ValueError("Missing required --output (peak-by-cell matrix path).")

ratio_dir = os.path.abspath(str(ratio_dir))
peakset   = os.path.abspath(str(peakset))
out_matrix= os.path.abspath(str(out_matrix))

if not os.path.isdir(ratio_dir):
    raise FileNotFoundError(f"--path not found or not a directory: {ratio_dir}")
if not os.path.exists(peakset):
    raise FileNotFoundError(f"--peakset not found: {peakset}")

WORKDIR = os.path.dirname(ratio_dir)
tmp_dir = os.path.join(WORKDIR, "tmp_peak2cell")

min_cover = int(config.get("min_cover", 3))

bedgraphs = sorted(glob.glob(os.path.join(ratio_dir, "*.bedgraph")))
if not bedgraphs:
    raise FileNotFoundError(f"No *.bedgraph found in: {ratio_dir}")

SAMPLES = [Path(x).name.replace(".bedgraph", "") for x in bedgraphs]


rule all:
    input:
        out_matrix

rule make_peak_id_bed:
    input:
        peakset
    output:
        peakid=temp(os.path.join(tmp_dir, "merge_peakID.DeChIC.bed"))
    shell:
        r"""
        mkdir -p $(dirname {output.peakid})
        awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,$1":"$2"-"$3}}' {input} > {output.peakid}
        """

rule peak2rate2cov_one_cell:
    input:
        peaks=os.path.join(tmp_dir, "merge_peakID.DeChIC.bed"),
        bg=lambda wc: os.path.join(ratio_dir, f"{wc.sample}.bedgraph")
    output:
        txt=temp(os.path.join(tmp_dir, "{sample}.peak2rate2cov.txt"))
    shell:
        r"""
        bedtools intersect -a {input.peaks} -b {input.bg} -wao \
          | awk 'BEGIN{{OFS="\t"}}{{print $4,$8,$9}}' \
          > {output.txt}
        """

rule build_peak_cell_matrix:
    input:
        peaks=os.path.join(tmp_dir, "merge_peakID.DeChIC.bed"),
        txts=expand(os.path.join(tmp_dir, "{sample}.peak2rate2cov.txt"), sample=SAMPLES)
    output:
        out=out_matrix
    params:
        min_cover=min_cover
    log:
        os.path.join(tmp_dir, "build_peak_cell_matrix.R.log")
    shell:
        r"""
set -euo pipefail
Rscript - <<'RS' > {log} 2>&1
suppressPackageStartupMessages(library(dplyr))

peaks_file <- "{input.peaks}"
out_file   <- "{output.out}"
min_cover  <- as.integer("{params.min_cover}")

peaks <- read.table(peaks_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
peak_order <- peaks$V4

files <- strsplit("{input.txts}", " +")[[1]]
files <- files[nchar(files) > 0]
sample_names <- sub("\\.peak2rate2cov\\.txt$", "", basename(files))

processed_list <- list()

for (k in seq_along(files)) {{
  f <- files[k]
  s <- sample_names[k]

  df <- read.table(f, sep="\t", header=FALSE, stringsAsFactors=FALSE,
                   na.strings=c("NA", ".", ""))
  colnames(df) <- c("peak","rate","covbp")

  df$rate  <- as.numeric(df$rate)
  df$covbp <- as.numeric(df$covbp)
  df$ind   <- ifelse(!is.na(df$covbp) & df$covbp > 0, 1, 0)

  tmp <- df %>%
    group_by(peak) %>%
    summarise(
      cover = sum(ind),
      avg_raw = ifelse(sum(ind)==0, NA_real_, mean(rate[ind==1], na.rm=TRUE))
    ) %>%
    mutate(
      avg = ifelse(cover < min_cover, NA_real_, avg_raw)
    ) %>%
    select(peak, avg)

  colnames(tmp) <- c("peak", s)
  processed_list[[s]] <- tmp
}}

mat <- Reduce(function(x,y) merge(x,y,by="peak",all=TRUE), processed_list)
mat <- mat[match(peak_order, mat$peak), ]

write.table(mat, file=out_file, sep="\t", quote=FALSE, row.names=FALSE)
RS
        """

