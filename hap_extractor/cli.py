"""Command-line interface definition."""

import argparse
import os

from . import __version__


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Generic GFF3 -> haplotype block extractor. "
            "If --all-del is omitted, a biallelic deletion TSV is first "
            "generated automatically from --vcf and then used in the analysis."
        )
    )

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    # Required inputs
    parser.add_argument(
        "--gff",
        required=True,
        help="Input GFF3 annotation file.",
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Input bgzip-compressed/indexed VCF file used for haplotype analysis.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help=(
            "Required output directory. All automatically named result files "
            "are written here; the directory is created automatically if it does not exist."
        ),
    )
    parser.add_argument(
        "--gene-list",
        required=True,
        help=(
            "Required gene selection. Use 'all' to process all genes in the GFF3 file, "
            "or provide a text file containing gene IDs, one ID per line."
        ),
    )

    # Optional all-DEL input. If absent, get_del() generates one from --vcf.
    parser.add_argument(
        "--all-del",
        default=None,
        help=(
            "Optional pre-generated deletion TSV with columns: "
            "chrom start end ref alt. If omitted, the script generates one "
            "from --vcf before the main analysis."
        ),
    )
    parser.add_argument(
        "--all-del-output",
        default=None,
        help=(
            "Output path for the automatically generated deletion TSV when "
            "--all-del is not supplied. Default: <output-dir>/all_dels.biallelic_only.<vcf-basename>.tsv"
        ),
    )
    parser.add_argument(
        "--report-every",
        type=int,
        default=100000,
        help="Progress interval while automatically generating all-DEL TSV (default: 100000).",
    )

    # Main analysis options
    parser.add_argument(
        "--mode",
        type=str.upper,
        choices=["GENE", "CDS"],
        default="CDS",
        help="Analysis mode: GENE or CDS (default: CDS).",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=2000,
        help="Flanking length for GENE mode (default: 2000).",
    )
    parser.add_argument(
        "--min-sample-count",
        type=int,
        default=2,
        help="Minimum haplotype occurrence count retained as a HAP row (default: 2).",
    )
    # Optional output filenames
    parser.add_argument(
        "--hap-output",
        default=None,
        help="Optional main output TSV filename/path. Default name is generated automatically.",
    )
    parser.add_argument(
        "--summary-output",
        default=None,
        help="Optional summary TSV filename/path. Default name is generated automatically.",
    )

    # Debug options
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Debug mode: only process the first --debug-n-genes genes.",
    )
    parser.add_argument(
        "--debug-n-genes",
        type=int,
        default=100,
        help="Number of genes processed in debug mode (default: 100).",
    )
    return parser


def parse_args(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.flank < 0:
        parser.error("--flank must be >= 0")
    if args.min_sample_count < 1:
        parser.error("--min-sample-count must be >= 1")
    if args.debug_n_genes < 1:
        parser.error("--debug-n-genes must be >= 1")
    if args.report_every < 1:
        parser.error("--report-every must be >= 1")

    if not os.path.isfile(args.gff):
        parser.error(f"GFF3 file not found: {args.gff}")
    if not os.path.isfile(args.vcf):
        parser.error(f"VCF file not found: {args.vcf}")
    if args.all_del and not os.path.isfile(args.all_del):
        parser.error(f"--all-del file not found: {args.all_del}")
    if args.gene_list.lower() != "all" and not os.path.isfile(args.gene_list):
        parser.error(
            "--gene-list must be 'all' or an existing gene-list file: "
            f"{args.gene_list}"
        )

    if args.gene_list.lower() != "all":
        has_gene_id = False
        with open(args.gene_list, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    has_gene_id = True
                    break
        if not has_gene_id:
            parser.error(f"--gene-list file contains no gene IDs: {args.gene_list}")

    return args