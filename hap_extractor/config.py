"""Runtime configuration object."""

from dataclasses import dataclass
import os
from typing import Optional


@dataclass(frozen=True)
class AnalysisConfig:
    gff: str
    vcf: str
    output_dir: str
    all_del: Optional[str] = None
    all_del_output: Optional[str] = None
    report_every: int = 100000
    mode: str = "CDS"
    flank: int = 2000
    min_sample_count: int = 1
    gene_list: str = "all"
    hap_output: Optional[str] = None
    summary_output: Optional[str] = None
    debug: bool = False
    debug_n_genes: int = 100

    @classmethod
    def from_namespace(cls, args):
        def expand(value):
            return os.path.expanduser(value) if value else None

        return cls(
            gff=expand(args.gff),
            vcf=expand(args.vcf),
            output_dir=expand(args.output_dir),
            all_del=expand(args.all_del),
            all_del_output=expand(args.all_del_output),
            report_every=args.report_every,
            mode=args.mode,
            flank=args.flank,
            min_sample_count=args.min_sample_count,
            gene_list=expand(args.gene_list),
            hap_output=args.hap_output,
            summary_output=args.summary_output,
            debug=args.debug,
            debug_n_genes=args.debug_n_genes,
        )