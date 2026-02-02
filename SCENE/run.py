# SCENE/run.py
import os
import argparse
import subprocess
from pathlib import Path
import textwrap

PKG_DIR = Path(__file__).resolve().parent  
MODULE_DIR = PKG_DIR / "modules"


class SCENEHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    """
    Combine:
      - ArgumentDefaultsHelpFormatter: show defaults in help (auto appends "(default: X)")
      - RawTextHelpFormatter: preserve manual newlines/indentation in help/epilog
    """
    pass


def _q(v):
    """Shell-safe-ish quoting for snakemake --config key=value ..."""
    v = str(v)
    if any(c in v for c in [' ', '\t', '\n', '"', "'"]):
        return '"' + v.replace('"', '\\"') + '"'
    return v


def run_snakemake(workflow: str, config_options: dict, cores: int = 4):
    snakefile = MODULE_DIR / f"{workflow}.snakefile"
    if not snakefile.exists():
        raise FileNotFoundError(f"Snakefile not found in package: {snakefile}")

    config_args = " ".join([f"{k}={_q(v)}" for k, v in config_options.items()])

    # Use absolute snakefile path so it works anywhere.
    default_config = PKG_DIR / "config.yaml"
    cmd = f"snakemake --snakefile {snakefile} --cores {cores} --configfile {default_config} --config {config_args}"
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def main():
    parser = argparse.ArgumentParser(
        description="SCENE is a single-cell DNA Deaminase based chromatin immuno-conversion sequencing data analysis pipeline",
        formatter_class=SCENEHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version", version="SCENE v1.0")

    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    # -----------------------------
    # rateMatrix
    # -----------------------------
    p_rate = subparsers.add_parser(
    "rateMatrix",
    help="Pre-processing to generate conversion-rate matrix of covered TC sites for each sample (and pseudo-bulk for single-cell mode).",
    formatter_class=SCENEHelpFormatter
)

    p_rate.add_argument(
        "-s", "--sample_type",
        choices=["bulk", "singlecell"],
        required=True,
        help="Sample type: bulk or singlecell"
    )
    p_rate.add_argument(
        "--path",
        required=True,
        metavar="DIR",
        help="Directory of raw fastq (*.R1.fastq.gz/*.R2.fastq.gz)"
    )
    p_rate.add_argument(
        "--SNP",
        required=True,
        metavar="BED",
        help="BED of CT/GA SNP sites (C->T and G->A)"
    )
    p_rate.add_argument(
        "--genome_fa",
        required=True,
        metavar="FASTA",
        help="Reference genome fasta used by BASAL (-d) and avgmod"
    )

    p_rate.add_argument(
        "--min_insert_size",
        type=int,
        default=1,
        metavar="<int>",
        help="-m  <int>    minimal insert size allowed"
    )
    p_rate.add_argument(
        "--max_insert_size",
        type=int,
        default=1000,
        metavar="<int>",
        help="-x  <int>    maximal insert size allowed"
    )
    p_rate.add_argument(
        "--report_repeat_hits",
        type=int,
        choices=[0, 1, 2],
        default=1,
        metavar="[0,1,2]",
        help=textwrap.dedent("""\
        -r  [0,1,2]  how to report repeat hits,
                    0=none(unique hit/pair); 1=random one; 2=all(slow)
        """).rstrip()
    )
    p_rate.add_argument(
        "--max_mismatch",
        type=float,
        default=0.06,
        metavar="<float>",
        help="-v  <float>  maximum percentage/number of mismatch bases in each read"
    )
    p_rate.add_argument(
        "--seed_size",
        type=int,
        default=16,
        metavar="<int>",
        help="-s  <int>    seed size (for WGBS/RRBS; min=8, max=16)"
    )
    p_rate.add_argument(
        "--random_seed",
        type=int,
        default=1,
        metavar="<int>",
        help="-S  <int>    seed for random number generation used in selecting multiple hits"
    )
    p_rate.add_argument(
        "--strand_info",
        type=int,
        choices=[0, 1, 2],
        default=1,
        metavar="[0,1,2]",
        help="-n  [0,1,2]  set mapping strand information"
    )
    p_rate.add_argument(
        "--max_gap",
        type=int,
        default=1,
        metavar="<int>",
        help="-g  <int>    maximum size of gap (deletion/insertion), <=3 bp"
    )
    p_rate.add_argument(
        "--convert",
        type=str,
        default="C:T",
        metavar="<str>",
        help="-M  <str>    convert-from and convert-to base(s) separated by ':', e.g. C:T"
    )

    # -----------------------------
    # callpeak
    # -----------------------------
    callpeak_epilog = textwrap.dedent("""\
    Notes:
      - The options starting with "met_" are passed to the metilene program.
      - filter/merge/blacklist parameters control peak filtering and merging.
    """)

    p_call = subparsers.add_parser(
        "callpeak",
        help="Call peaks for each sample (and pseudo-bulk in single-cell mode).",
        formatter_class=SCENEHelpFormatter,
        epilog=callpeak_epilog
    )
    p_call.add_argument(
        "-s", "--sample_type",
        choices=["bulk", "singlecell"],
        required=True,
        help="Sample type: bulk or singlecell"
    )
    p_call.add_argument(
        "--path",
        required=True,
        metavar="DIR",
        help="Directory of bedgraphs (ratio/)."
    )
    p_call.add_argument(
        "-IgG",
        dest="IgG",
        default=None,
        metavar="BEDGRAPH",
        help="Optional: shared IgG bedgraph path. If omitted, use mean/pseudoIgG mode."
    )

    p_call.add_argument(
        "--met_M",
        type=int,
        default=50,
        metavar="<n>",
        help="-M, --maxdist <n>      maximum distance"
    )
    p_call.add_argument(
        "--met_m",
        type=int,
        default=6,
        metavar="<n>",
        help="-m, --mincpgs <n>      minimum CpGs"
    )
    p_call.add_argument(
        "--met_t",
        type=int,
        default=40,
        metavar="<n>",
        help="-t, --threads <n>      number of threads"
    )

    # ---- filtering/merging/blacklist params
    p_call.add_argument(
        "--p_th",
        type=float,
        default=0.05,
        metavar="<float>",
        help="Maximum p-value (<) for output peaks"
    )
    p_call.add_argument(
        "--d_th",
        type=float,
        default=0.2,
        metavar="<float>",
        help="Minimum mean methylation difference (>=)"
    )
    p_call.add_argument(
        "--merge_d",
        type=int,
        default=200,
        metavar="<int>",
        help="Minimum distance for merging nearby peaks"
    )
    p_call.add_argument(
        "--bl_f",
        type=float,
        default=0.5,
        metavar="<float>",
        help="Minimum overlap fraction with blacklist; peaks exceeding this fraction will be removed"
    )

    # -----------------------------
    # cellMatrix
    # -----------------------------
    p_p2c = subparsers.add_parser(
        "cellMatrix",
        help="Build a peak-by-cell matrix from the pseudo-bulk peak set.",
        formatter_class=SCENEHelpFormatter
    )
    p_p2c.add_argument("--path", required=True, metavar="DIR", help="Directory of per-cell bedgraphs (ratio/).")
    p_p2c.add_argument("--peakset", required=True, metavar="BED", help="Merged peak bed file (e.g. merged callpeak output).")
    p_p2c.add_argument("--output", required=True, metavar="FILE", help="Output peak-by-cell matrix file (txt).")

    # -----------------------------
    # peakPresence
    # -----------------------------
    p_pp = subparsers.add_parser(
        "peakPresence",
        help="Compute Peak Presence Percentage.",
        formatter_class=SCENEHelpFormatter
    )
    p_pp.add_argument("--path", required=True, metavar="FILE", help="Peak-by-cell matrix file (txt).")
    p_pp.add_argument("--output", required=True, metavar="FILE", help="Output peakPresence matrix file (txt).")

    args = parser.parse_args()

    if args.command == "rateMatrix":
        config_options = {
            "sample_type": args.sample_type,
            "raw_dir": os.path.abspath(args.path),
            "snp_bed": os.path.abspath(args.SNP),
            "genome_fa": os.path.abspath(args.genome_fa),
            "align_m": args.min_insert_size,
            "align_x": args.max_insert_size,
            "align_r": args.report_repeat_hits,
            "align_v": args.max_mismatch,
            "align_s": args.seed_size,
            "align_S": args.random_seed,
            "align_n": args.strand_info,
            "align_g": args.max_gap,
            "align_M": args.convert,
        }
        run_snakemake("rateMatrix", config_options)

    elif args.command == "callpeak":
        config_options = {
            "sample_type": args.sample_type,
            "path": os.path.abspath(args.path),
        }
        if args.IgG is not None:
            config_options["IgG_bedgraph"] = os.path.abspath(args.IgG)

        config_options.update({
            "met_M": args.met_M,
            "met_m": args.met_m,
            "met_t": args.met_t,
            "p_th": args.p_th,
            "d_th": args.d_th,
            "merge_d": args.merge_d,
            "bl_f": args.bl_f,
        })

        run_snakemake("callpeak", config_options)

    elif args.command == "cellMatrix":
        config_options = {
            "path": os.path.abspath(args.path),
            "peakset": os.path.abspath(args.peakset),
            "output": os.path.abspath(args.output),
        }
        run_snakemake("cellMatrix", config_options)

    elif args.command == "peakPresence":
        config_options = {
            "path": os.path.abspath(args.path),
            "output": os.path.abspath(args.output),
        }
        run_snakemake("peakPresence", config_options)

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
