# modules/peakPresence.snakefile
from pathlib import Path
import os

BASE = Path(workflow.basedir)
PKG  = BASE.parent
configfile: str(PKG / "config.yaml")

peak2cell = config.get("path", None)      # --path
out_matrix = config.get("output", None)   # --output

if peak2cell is None or str(peak2cell).strip() == "":
    raise ValueError("Missing required --path (peak-by-cell matrix).")
if out_matrix is None or str(out_matrix).strip() == "":
    raise ValueError("Missing required --output (peakPresenceMatrix.txt).")

peak2cell = os.path.abspath(str(peak2cell))
out_matrix = os.path.abspath(str(out_matrix))

if not os.path.exists(peak2cell):
    raise FileNotFoundError(f"--path not found: {peak2cell}")

min_detect_prop  = float(config.get("min_detect_prop", 0.5))
rate_cutoff      = float(config.get("rate_cutoff", 0.15))
max_tcga_density = float(config.get("max_tcga_density", 0.2))


rule all:
    input:
        out_matrix,

rule peak_presence:
    input:
        peak2cell
    output:
        matrix=out_matrix,
    log:
        "logs/peakPresence.R.log"
    shell:
        r"""
        mkdir -p logs
        /mnt/datadisk/chenxiyang/miniconda3/envs/NChIP/bin/Rscript - <<'RS' > {log} 2>&1
        suppressPackageStartupMessages({{
          library(BSgenome)
          library(Biostrings)
          library(GenomicRanges)
          library(BSgenome.Mmusculus.UCSC.mm10)
        }})

        in_file <- "{input}"
        out_mat <- "{output.matrix}"

        min_detect_prop  <- as.numeric("{min_detect_prop}")
        rate_cutoff      <- as.numeric("{rate_cutoff}")
        max_tcga_density <- as.numeric("{max_tcga_density}")

        peak_cell_matrix <- read.table(in_file, sep="\t", header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)

        non_na_prop <- rowMeans(!is.na(peak_cell_matrix))
        mat1 <- peak_cell_matrix[non_na_prop >= min_detect_prop, , drop=FALSE]
        if (nrow(mat1) == 0) stop("No peaks remain after non-NA proportion filter.")

        parse_peak <- function(x){{
          parts1 <- strsplit(x, ":", fixed=TRUE)
          chr <- vapply(parts1, function(z) z[1], character(1))
          se  <- vapply(parts1, function(z) z[2], character(1))
          parts2 <- strsplit(se, "-", fixed=TRUE)
          start <- as.integer(vapply(parts2, function(z) z[1], character(1)))
          end   <- as.integer(vapply(parts2, function(z) z[2], character(1)))
          data.frame(chr=chr, start=start, end=end, stringsAsFactors=FALSE)
        }}

        df_peaks <- parse_peak(rownames(mat1))
        peaks_gr <- GRanges(seqnames=df_peaks$chr, ranges=IRanges(start=df_peaks$start, end=df_peaks$end))

        seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, peaks_gr)
        tc_counts <- vcountPattern("TC", seqs)
        ga_counts <- vcountPattern("GA", seqs)
        peak_lengths <- width(seqs)
        tcga_density <- (tc_counts + ga_counts) / peak_lengths

        keep2 <- tcga_density < max_tcga_density
        mat2 <- mat1[keep2, , drop=FALSE]
        if (nrow(mat2) == 0) stop("No peaks remain after TC+GA density filter (< threshold).")

        bin <- as.data.frame(lapply(mat2, function(x){{
          ifelse(is.na(x), NA_real_, ifelse(x > rate_cutoff, 1, 0))
        }}), check.names=FALSE)
        rownames(bin) <- rownames(mat2)

        calc_prop <- function(x){{
          v <- x[!is.na(x)]
          sum(v == 1) / length(v)
        }}
        prop <- apply(bin, 1, calc_prop)

        df2 <- parse_peak(names(prop))
        out_df <- data.frame(chr=df2$chr, start=df2$start, end=df2$end, PeakPresence=as.numeric(prop), stringsAsFactors=FALSE)
        write.table(out_df, file=out_mat, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
RS
        """
