# modules/callpeak.snakefile
from pathlib import Path
import os, glob

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


# -----------------------------
# required CLI-injected config
# -----------------------------
sample_type = config.get("sample_type", None)
if sample_type not in ("bulk", "singlecell"):
    raise ValueError("Missing/invalid sample_type. Use -s bulk or -s singlecell.")

ratio_dir = config.get("path", None)  # --path
if ratio_dir is None or str(ratio_dir).strip() == "":
    raise ValueError("Missing required --path (SCENE callpeak --path BEDGRAPH_DIR ...)")
ratio_dir = os.path.abspath(str(ratio_dir))
if not os.path.isdir(ratio_dir):
    raise FileNotFoundError(f"--path not found or not a directory: {ratio_dir}")

WORKDIR = os.path.dirname(ratio_dir)

# output dirs (project-relative)
call_dir = os.path.join(WORKDIR, config.get("paths", {}).get("callpeak_dir", "callpeak"))
merged_ratio_dir = os.path.join(WORKDIR, config.get("paths", {}).get("merged_ratio_dir", "merged_bam_ratio"))

# package resources/scripts
blacklist_bed = pkg_file(config.get("paths", {}).get("blacklist_bed", "resources/mm10.DeChIC.blacklist.bed"))
metilene_input_pl = pkg_file(config.get("paths", {}).get("metilene_input_pl", "scripts/metilene_input.pl"))

for fp, label in [
    (blacklist_bed, "blacklist_bed"),
    (metilene_input_pl, "metilene_input_pl"),
]:
    if fp == "" or (not Path(fp).exists()):
        raise FileNotFoundError(f"{label} not found in package: {fp}")

# IgG optional
igg_bedgraph = config.get("IgG_bedgraph", None)
use_igg = (igg_bedgraph is not None) and (str(igg_bedgraph).strip() != "")
igg_path = os.path.abspath(str(igg_bedgraph)) if use_igg else ""
if use_igg and (not Path(igg_path).exists()):
    raise FileNotFoundError(f"IgG bedgraph not found: {igg_path}")


# -----------------------------
# metilene/filter params (from config.yaml)
# -----------------------------
met_M = int(config.get("metilene", {}).get("M", 50))
met_m = int(config.get("metilene", {}).get("m", 6))
met_t = int(config.get("metilene", {}).get("t", 40))

p_th    = float(config.get("filter", {}).get("p_threshold", 0.05))
d_th    = float(config.get("filter", {}).get("d_threshold", 0.2))
merge_d = int(config.get("filter", {}).get("merge_distance", 200))
bl_f    = float(config.get("filter", {}).get("blacklist_fraction", 0.5))

# -----------------------------
# discover samples
# -----------------------------
indiv_bedgraphs = sorted(glob.glob(os.path.join(ratio_dir, "*.bedgraph")))
if not indiv_bedgraphs:
    raise FileNotFoundError(f"No *.bedgraph found in --path: {ratio_dir}")

INDIV_SAMPLES = [Path(x).name.replace(".bedgraph", "") for x in indiv_bedgraphs]

MERGED_SAMPLES = []
if sample_type == "singlecell":
    merged_bedgraphs = sorted(glob.glob(os.path.join(merged_ratio_dir, "*.merge.bedgraph")))
    if not merged_bedgraphs:
        raise FileNotFoundError(
            f"sample_type=singlecell but no merged *.merge.bedgraph found in {merged_ratio_dir}. "
            f"Please run: SCENE rateMatrix -s singlecell ... first (merge is generated there)."
        )
    MERGED_SAMPLES = [Path(x).name.replace(".bedgraph", "") for x in merged_bedgraphs]

SAMPLES = sorted(list(dict.fromkeys(INDIV_SAMPLES + MERGED_SAMPLES)))

def bedgraph_path(wc):
    if str(wc.sample).endswith(".merge"):
        return os.path.join(merged_ratio_dir, f"{wc.sample}.bedgraph")
    return os.path.join(ratio_dir, f"{wc.sample}.bedgraph")

# -----------------------------
# targets
# -----------------------------
rule all:
    input:
        expand(os.path.join(call_dir, "{sample}.DeChIC.bed"), sample=SAMPLES)

# -----------------------------
# pseudo IgG
# -----------------------------
rule make_pseudo_igg:
    input:
        bg=bedgraph_path
    output:
        pseudo=temp(os.path.join(call_dir, "{sample}.pseudoIgG.bedgraph"))
    shell:
        r"""
        mkdir -p $(dirname {output.pseudo})
        awk 'BEGIN{{OFS="\t"}}
             NR==FNR{{sum+=$4; n+=1; next}}
             FNR==1{{mean=sum/n}}
             {{print $1,$2,$3,mean}}
            ' {input.bg} {input.bg} > {output.pseudo}
        """

# -----------------------------
# metilene input
# -----------------------------
rule metilene_input:
    input:
        in1=bedgraph_path,
        in2=lambda wc: igg_path if use_igg else os.path.join(call_dir, f"{wc.sample}.pseudoIgG.bedgraph")
    output:
        out=os.path.join(call_dir, "{sample}.input.txt")
    shell:
        r"""
        mkdir -p $(dirname {output.out})
        perl {metilene_input_pl} \
          --in1 {input.in1} \
          --in2 {input.in2} \
          --h1 target \
          --h2 control \
          --out {output.out}
        """

# -----------------------------
# metilene
# -----------------------------
rule metilene_one:
    input:
        txt=os.path.join(call_dir, "{sample}.input.txt")
    output:
        out=os.path.join(call_dir, "{sample}.output.txt")
    log:
        os.path.join(call_dir, "logs", "{sample}.metilene.log")
    shell:
        r"""
        mkdir -p $(dirname {log})
        metilene -a "target" -b "control" -M {met_M} -m {met_m} -t {met_t} \
          {input.txt} > {output.out} 2> {log}
        """

# -----------------------------
# filter -> DeChIC.bed
# -----------------------------
rule callpeak_filter:
    input:
        out=os.path.join(call_dir, "{sample}.output.txt")
    output:
        bed=os.path.join(call_dir, "{sample}.DeChIC.bed")
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        awk -v p="{p_th}" -v d="{d_th}" '($7+0 < p && $5+0 > d){{print}}' {input.out} \
          | cut -f 1,2,3 \
          | sort -k1,1 -k2,2n \
          | bedtools merge -i - -d {merge_d} \
          | bedtools intersect -v -a - -b {blacklist_bed} -f {bl_f} \
          | uniq > {output.bed}
        """

