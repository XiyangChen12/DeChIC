#!/usr/bin/python3
# -*- coding:utf-8 -*-

import glob
import os
import re
from pathlib import Path

# ============================================================
# Package paths
# ============================================================
BASE = Path(workflow.basedir)
PKG  = BASE.parent 

configfile: str(PKG / "config.yaml")


def pkg_file(p: str) -> str:
    if p is None:
        return ""
    p = str(p).strip()
    if p == "":
        return ""
    pp = Path(p)
    return str(pp if pp.is_absolute() else (PKG / pp))


def maybe_cmd_or_path(p: str) -> str:
    if p is None:
        return ""
    p = str(p).strip()
    if p == "":
        return ""
    if ("/" in p) or ("\\" in p) or p.startswith("."):
        return pkg_file(p)
    return p


# -----------------------------
# required CLI-injected config
# -----------------------------
sample_type = config.get("sample_type", None)
if sample_type not in ("bulk", "singlecell"):
    raise ValueError("Missing/invalid sample_type. Use -s bulk or -s singlecell.")

raw_dir = config.get("raw_dir", None)
if raw_dir is None or str(raw_dir).strip() == "":
    raise ValueError("Missing required config key: raw_dir (SCENE rateMatrix --path FASTQ_DIR ...)")
raw_dir = os.path.abspath(str(raw_dir))
if not os.path.isdir(raw_dir):
    raise FileNotFoundError(f"raw_dir not found or not a directory: {raw_dir}")

snp_bed = config.get("snp_bed", None)
if snp_bed is None or str(snp_bed).strip() == "":
    raise ValueError("Missing required config key: snp_bed (SCENE rateMatrix --SNP CT_GA.snp.bed ...)")
snp_bed = os.path.abspath(str(snp_bed))
if not os.path.exists(snp_bed):
    raise FileNotFoundError(f"snp_bed not found: {snp_bed}")

genome_fa = config.get("genome_fa", None)
if genome_fa is None or str(genome_fa).strip() == "":
    raise ValueError("Missing required config key: genome_fa (SCENE rateMatrix --genome_fa mm10.fa ...)")
genome_fa = os.path.abspath(str(genome_fa))
if not os.path.exists(genome_fa):
    raise FileNotFoundError(f"genome_fa not found: {genome_fa}")

# -----------------------------
# sample_type specific params
# -----------------------------
if sample_type == "bulk":
    clip_R1 = config["bulk"]["clip_R1"]
    clip_R2 = config["bulk"]["clip_R2"]
    three_prime_clip_R1 = config["bulk"]["three_prime_clip_R1"]
    three_prime_clip_R2 = config["bulk"]["three_prime_clip_R2"]
    depth = config["bulk"]["depth"]
else:
    clip_R1 = config["singlecell"]["clip_R1"]
    clip_R2 = config["singlecell"]["clip_R2"]
    three_prime_clip_R1 = config["singlecell"]["three_prime_clip_R1"]
    three_prime_clip_R2 = config["singlecell"]["three_prime_clip_R2"]
    depth = config["singlecell"]["depth"]

# -----------------------------
# tools / package resources
# -----------------------------
basal_bin = maybe_cmd_or_path(config.get("paths", {}).get("basal_bin", "basal"))
basalkit_avgmod_py = pkg_file(config.get("paths", {}).get("basalkit_avgmod_py", "scripts/basalkit_filterTC.py"))
blacklist_bed = pkg_file(config.get("paths", {}).get("blacklist_bed", "resources/mm10.DeChIC.blacklist.bed"))

if basalkit_avgmod_py == "" or not Path(basalkit_avgmod_py).exists():
    raise FileNotFoundError(f"basalkit_avgmod_py not found in package: {basalkit_avgmod_py}")
if blacklist_bed == "" or not Path(blacklist_bed).exists():
    raise FileNotFoundError(f"blacklist_bed not found in package: {blacklist_bed}")

# output dirs (project-relative, DO NOT resolve to package)
merged_bam_dir = config.get("paths", {}).get("merged_bam_dir", "merged_bam_file")
merged_ratio_dir = config.get("paths", {}).get("merged_ratio_dir", "merged_bam_ratio")

# -----------------------------
# BASAL align params (CLI overrides if injected; else config.yaml defaults)
# -----------------------------
def _align_param(flat_key, yaml_key):
    if flat_key in config and config[flat_key] is not None and str(config[flat_key]).strip() != "":
        return config[flat_key]
    return config["basal"]["align"][yaml_key]

align_m = _align_param("align_m", "min_insert_size")
align_x = _align_param("align_x", "max_insert_size")
align_r = _align_param("align_r", "report_repeat_hits")
align_v = _align_param("align_v", "max_mismatch")
align_s = _align_param("align_s", "seed_size")
align_S = _align_param("align_S", "random_seed")
align_n = _align_param("align_n", "strand_info")
align_g = _align_param("align_g", "max_gap")
align_M = _align_param("align_M", "convert")

# -----------------------------
# discover samples
# -----------------------------
SAMPLES = sorted([
    os.path.basename(f).replace(".R1.fastq.gz", "")
    for f in glob.glob(os.path.join(raw_dir, "*.R1.fastq.gz"))
])
if not SAMPLES:
    raise FileNotFoundError(f"No *.R1.fastq.gz found in raw_dir: {raw_dir}")

# singlecell groups: <GROUP>_<digits>
GROUPS = []
if sample_type == "singlecell":
    pat = re.compile(r"^(.*)_(\d+)$")
    groups = []
    for sname in SAMPLES:
        m = pat.match(sname)
        if m:
            groups.append(m.group(1))
    GROUPS = sorted(list(dict.fromkeys(groups)))
    if not GROUPS:
        raise ValueError("sample_type=singlecell but no samples match '<GROUP>_<digits>' from FASTQ names.")

print("raw_dir:", raw_dir)
print("sample_type:", sample_type)
print("SAMPLES:", SAMPLES)
print("GROUPS:", GROUPS if sample_type == "singlecell" else "NA")
print("snp_bed:", snp_bed)
print("genome_fa:", genome_fa)
print("basal_bin:", basal_bin)
print("basalkit_avgmod_py:", basalkit_avgmod_py)
print("blacklist_bed:", blacklist_bed)

# -----------------------------
# rule all
# -----------------------------
rule all:
    input:
        expand("ratio/{sample}.bedgraph", sample=SAMPLES),
        expand(os.path.join(merged_ratio_dir, "{group}.merge.bedgraph"), group=GROUPS) if sample_type == "singlecell" else []

# -----------------------------
# trim_galore
# -----------------------------
rule filter:
    input:
        R1=lambda wc: os.path.join(raw_dir, f"{wc.sample}.R1.fastq.gz"),
        R2=lambda wc: os.path.join(raw_dir, f"{wc.sample}.R2.fastq.gz")
    output:
        R1="{sample}.R1_val_1.fq.gz",
        R2="{sample}.R2_val_2.fq.gz"
    log:
        "logs/trim_galore/{sample}.log"
    threads: 7
    shell:
        r"""
        mkdir -p logs/trim_galore
        trim_galore --trim-n \
          --clip_R1 {clip_R1} --clip_R2 {clip_R2} \
          --three_prime_clip_R1 {three_prime_clip_R1} \
          --three_prime_clip_R2 {three_prime_clip_R2} \
          --length 50 -q 20 --paired -j {threads} \
          {input.R1} {input.R2} -o ./ > {log} 2>&1
        """

# -----------------------------
# align (per-cell)
# -----------------------------
rule align:
    input:
        R1="{sample}.R1_val_1.fq.gz",
        R2="{sample}.R2_val_2.fq.gz"
    output:
        "bam_file/{sample}.bam"
    log:
        "logs/basal/{sample}.log"
    threads: 64
    params:
        m=align_m,
        x=align_x,
        r=align_r,
        v=align_v,
        s=align_s,
        S=align_S,
        n=align_n,
        g=align_g,
        M=align_M
    shell:
        r"""
        mkdir -p bam_file logs/basal
        {basal_bin} -a {input.R1} -b {input.R2} -d {genome_fa} \
          -m {params.m} -x {params.x} -p {threads} -r {params.r} \
          -v {params.v} -s {params.s} -S {params.S} \
          -n {params.n} -g {params.g} -M {params.M} \
          -o {output} 2> {log}
        """

# -----------------------------
# sort bam (per-cell)
# -----------------------------
rule sortBam:
    input:
        "bam_file/{sample}.bam"
    output:
        "bam_file/{sample}.sort.bam"
    threads: 64
    shell:
        r"""
        mkdir -p bam_file
        samtools sort -@ {threads} -o {output} {input}
        """

# -----------------------------
# basalkit avgmod (per-cell)
# -----------------------------
rule ratio:
    input:
        "bam_file/{sample}.sort.bam"
    output:
        "{sample}_AvgMod.tsv"
    shell:
        r"""
        python3 {basalkit_avgmod_py} avgmod {input} {genome_fa} \
          -m 1 -o {wildcards.sample}

        # normalize
        if [ -f "{wildcards.sample}_AvgMod.tsv" ] && [ ! -f "{output}" ]; then
          mv "{wildcards.sample}_AvgMod.tsv" "{output}"
        fi
        """

# -----------------------------
# capture TC (per-cell)
# -----------------------------
rule captureTC:
    input:
        "{sample}_AvgMod.tsv"
    output:
        temp("ratio/{sample}.TC.txt")
    shell:
        r"""
        mkdir -p ratio
        awk 'BEGIN{{OFS="\t"}} NR>1 && $4=="TC" && $6 + 0 >{depth} && $5!="NA" \
        {{print $1,$2-1,$2,1-$5}}' {input} > {output}
        """

# -----------------------------
# remove SNP (CT/GA) + blacklist -> bedgraph (per-cell)
# -----------------------------
rule removeSNP_blacklist:
    input:
        "ratio/{sample}.TC.txt"
    output:
        "ratio/{sample}.bedgraph"
    params:
        SNP=snp_bed,
        BL=blacklist_bed
    shell:
        r"""
        mkdir -p ratio
        bedtools intersect -v -a {input} -b {params.SNP} \
        | bedtools intersect -v -a - -b {params.BL} \
        | cut -f 1-4 \
        | sort -k1,1 -k2,2n > {output}
        """

# ============================================================
# pseudo-bulk if in singlecell mode
# ============================================================
if sample_type == "singlecell":

    def group_bams(wc):
        pat = re.compile(rf"^{re.escape(wc.group)}_(\d+)\.bam$")
        bams = []
        for sname in SAMPLES:
            if pat.match(f"{sname}.bam"):
                bams.append(f"bam_file/{sname}.bam")
        if not bams:
            bams = sorted(glob.glob(f"bam_file/{wc.group}_[0-9]*.bam"))
        if not bams:
            raise ValueError(f"No BAMs found for group={wc.group} under bam_file/. Check FASTQ naming '<GROUP>_<digits>'.")
        return bams

    rule merge_group_bam:
        input:
            group_bams
        output:
            os.path.join(merged_bam_dir, "{group}.merge.bam")
        threads: 16
        shell:
            r"""
            mkdir -p {merged_bam_dir}
            samtools merge -@ {threads} {output} {input}
            """

    rule sort_merged_bam:
        input:
            os.path.join(merged_bam_dir, "{group}.merge.bam")
        output:
            os.path.join(merged_bam_dir, "{group}.merge.sort.bam")
        threads: 16
        shell:
            r"""
            mkdir -p {merged_bam_dir}
            samtools sort -@ {threads} -o {output} {input}
            """

    rule ratio_merged:
        input:
            os.path.join(merged_bam_dir, "{group}.merge.sort.bam")
        output:
            os.path.join(merged_ratio_dir, "{group}.merge_AvgMod.tsv")
        shell:
            r"""
            mkdir -p {merged_ratio_dir}
            python3 {basalkit_avgmod_py} avgmod {input} {genome_fa} -m 2 -o {wildcards.group}.merge

            # normalize
            if [ -f "{wildcards.group}.merge_AvgMod.tsv" ] && [ ! -f "{output}" ]; then
              mv "{wildcards.group}.merge_AvgMod.tsv" "{output}"
            fi
            """

    rule captureTC_merged:
        input:
            "merged_bam_ratio/{group}.merge_AvgMod.tsv"
        output:
            temp("{merged_ratio_dir}/{group}.merge.TC.txt")
        shell:
            r"""
            mkdir -p {merged_ratio_dir}
            awk 'BEGIN{{OFS="\t"}} NR>1 && $4=="TC" && $5!="NA" \
            {{print $1,$2-1,$2,1-$5}}' {input} > {output}
            """

    rule removeSNP_blacklist_merged:
        input:
            os.path.join(merged_ratio_dir, "{group}.merge.TC.txt")
        output:
            os.path.join(merged_ratio_dir, "{group}.merge.bedgraph")
        params:
            SNP=snp_bed,
            BL=blacklist_bed
        shell:
            r"""
            mkdir -p {merged_ratio_dir}
            bedtools intersect -v -a {input} -b {params.SNP} \
            | bedtools intersect -v -a - -b {params.BL} \
            | cut -f 1-4 \
            | sort -k1,1 -k2,2n > {output}
            """

