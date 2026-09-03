"""Main analysis pipeline."""

from collections import defaultdict
import os

import pysam

from .constants import META_QUERY_REGION
from .deletions import (
    extract_cds_overlapping_dels_from_list,
    extract_gene_overlapping_dels_from_list,
    get_del,
    load_all_dels_tsv,
)
from .gene_filter import build_gene_id_filter, gene_id_in_list
from .gff3 import load_gene_regions
from .haplotypes import build_haplotype_rows, cluster_all_haplotypes
from .large_output import compress_large_haplotype_tsv_if_needed
from .output import (
    make_cds_ranges_array,
    write_block_and_hap_rows,
    write_format_header,
    write_meta_rows,
)
from .paths import (
    ensure_parent_dir,
    load_gene_ids_from_list,
    make_default_all_del_path,
    make_default_output_names,
    make_unique_path,
    resolve_output_path,
)
from .vcf_utils import (
    extract_cds_variants,
    extract_region_variants,
    fetch_site_genotypes_from_vcf,
)


def run_pipeline(config):
    """Run the complete GFF3/VCF -> haplotype-block workflow."""
    output_dir = config.output_dir
    os.makedirs(output_dir, exist_ok=True)

    gene_id_list = load_gene_ids_from_list(config.gene_list)

    default_output_file, default_summary_file = make_default_output_names(
        config.mode,
        config.debug,
        config.min_sample_count,
        config.vcf,
    )

    output_path = resolve_output_path(
        output_dir,
        config.hap_output,
        default_output_file,
    )
    output_path = make_unique_path(
        output_path,
        related_suffixes=(".gz", ".gz.tbi"),
    )

    summary_path = resolve_output_path(
        output_dir,
        config.summary_output,
        default_summary_file,
    )
    summary_path = make_unique_path(summary_path)

    ensure_parent_dir(output_path)
    ensure_parent_dir(summary_path)

    # --all-del is optional:
    #   supplied     -> use it directly
    #   not supplied -> generate it first from --vcf using get_del()
    if config.all_del:
        del_position_file = config.all_del
        print(f"[INFO] Using provided all-DEL TSV: {del_position_file}")
    else:
        if config.all_del_output:
            del_position_file = config.all_del_output
        else:
            del_position_file = make_default_all_del_path(
                config.vcf,
                output_dir,
            )

        del_position_file = make_unique_path(del_position_file)

        get_del(
            config.vcf,
            del_position_file,
            report_every=config.report_every,
        )

    print("Loading deletion TSV...")
    dels_by_chr = load_all_dels_tsv(del_position_file)

    print("Loading VCF...")
    vcf = pysam.VariantFile(config.vcf)
    samples = list(vcf.header.samples)
    total_samples = len(samples)

    print("Building gene ID filter...")
    raw_gene_ids, match_gene_ids = build_gene_id_filter(gene_id_list)

    if match_gene_ids is None:
        print("Gene selection: all genes (--gene-list all)")
    else:
        print(
            "Gene list provided:",
            len(raw_gene_ids),
            "raw IDs",
        )

    print("Loading and normalizing GFF3...")
    (
        genes,
        transcript_units_map,
        exon_units_map,
    ) = load_gene_regions(config.gff)

    total_genes_in_gff = len(genes)

    if match_gene_ids is not None:
        genes = [
            g
            for g in genes
            if gene_id_in_list(
                g[0],
                match_gene_ids,
            )
        ]

        print(
            "Gene ID filtering:",
            total_genes_in_gff,
            "->",
            len(genes),
            "matched genes",
        )

    if config.debug:
        genes = genes[:config.debug_n_genes]

    gene_count = len(genes)

    print("VCF loaded with", total_samples, "samples")
    print("Total genes to process:", gene_count)
    print("MODE =", config.mode)
    print(
        "Transcript behavior in CDS mode: "
        "ALL normalized transcripts with CDS"
    )

    processed = 0
    total_hap_count = 0
    total_blocks_written = 0
    total_transcript_blocks_written = 0

    total_transcripts_seen = 0
    total_transcripts_with_cds = 0
    skipped_transcript_no_cds = 0
    hap_sample_count_map = defaultdict(int)

    def finalize_and_write_block(
        out,
        block_chrom,
        block_start_1based,
        block_end_1based,
        block_id,
        gene_id,
        strand,
        block_mode,
        flank_value,
        cds_ranges_array,
        region,
        transcript_id=".",
    ):
        nonlocal total_hap_count
        nonlocal total_blocks_written
        nonlocal total_transcript_blocks_written

        hap_dict = cluster_all_haplotypes(
            region,
            samples,
        )

        (
            hap_rows,
            kept_hap_count,
            hap_updates,
        ) = build_haplotype_rows(
            hap_dict,
            config.min_sample_count,
        )

        total_hap_count += kept_hap_count
        total_blocks_written += 1

        if transcript_id != ".":
            total_transcript_blocks_written += 1

        for sc, n in hap_updates.items():
            hap_sample_count_map[sc] += n

        write_block_and_hap_rows(
            out=out,
            block_chrom=block_chrom,
            block_start_1based=block_start_1based,
            block_end_1based=block_end_1based,
            block_id=block_id,
            gene_id=gene_id,
            strand=strand,
            block_mode=block_mode,
            flank_value=flank_value,
            cds_ranges_array=cds_ranges_array,
            region=region,
            hap_rows=hap_rows,
        )

    try:
        with open(
            output_path,
            "w",
            encoding="utf-8",
        ) as out:
            write_format_header(out)
            write_meta_rows(
                out,
                total_samples,
                samples,
                config.min_sample_count,
                config.mode,
                gene_id_list,
            )

            for (
                gene_id,
                chrom,
                gstart,
                gend,
                strand,
            ) in genes:
                processed += 1
                print(
                    f"Processing {processed}/{gene_count}: "
                    f"{gene_id}"
                )

                # ------------------------------------------------
                # GENE MODE
                # ------------------------------------------------
                if config.mode == "GENE":
                    region_start = max(
                        1,
                        gstart - config.flank,
                    )
                    region_end = gend + config.flank

                    region = extract_region_variants(
                        vcf,
                        chrom,
                        region_start - 1,
                        region_end,
                    )

                    dels_chr = dels_by_chr.get(
                        chrom,
                        [],
                    )

                    extra_dels = (
                        extract_gene_overlapping_dels_from_list(
                            dels_chr,
                            region_start,
                            region_end,
                        )
                    )

                    var_map = {
                        (
                            s["pos"],
                            s["ref"],
                            s["alt"],
                        ): s
                        for s in region
                    }

                    for s in extra_dels:
                        k = (
                            s["pos"],
                            s["ref"],
                            s["alt"],
                        )

                        if k in var_map:
                            continue

                        geno = fetch_site_genotypes_from_vcf(
                            vcf,
                            chrom,
                            s["pos"],
                            s["ref"],
                            s["alt"],
                        )

                        if geno is None:
                            continue

                        s2 = dict(s)
                        s2["genotypes"] = geno
                        var_map[k] = s2

                    region = [
                        var_map[k]
                        for k in sorted(
                            var_map,
                            key=lambda x: x[0],
                        )
                    ]

                    if not region:
                        print(
                            f"[SKIP] {gene_id}: "
                            f"no variants found in "
                            f"{chrom}:{region_start}-{region_end}"
                        )
                        continue

                    block_chrom = chrom
                    block_start_1based = region_start
                    block_end_1based = region_end
                    block_id = (
                        f"{block_chrom}:"
                        f"{block_start_1based}-"
                        f"{block_end_1based}|"
                        f"{gene_id}"
                    )

                    finalize_and_write_block(
                        out=out,
                        block_chrom=block_chrom,
                        block_start_1based=block_start_1based,
                        block_end_1based=block_end_1based,
                        block_id=block_id,
                        gene_id=gene_id,
                        strand=strand,
                        block_mode="GENE",
                        flank_value=config.flank,
                        cds_ranges_array=[],
                        region=region,
                    )

                # ------------------------------------------------
                # CDS MODE: ALL TRANSCRIPTS
                # ------------------------------------------------
                elif config.mode == "CDS":
                    tx_units = transcript_units_map.get(
                        gene_id,
                        [],
                    )

                    if not tx_units:
                        print(
                            f"[SKIP] {gene_id}: "
                            "no normalized transcript units found"
                        )
                        continue

                    total_transcripts_seen += len(tx_units)

                    for tx_i, tx_unit in enumerate(
                        tx_units,
                        start=1,
                    ):
                        tx_id = (
                            tx_unit.get("transcript_id")
                            or "."
                        )
                        cds_ranges = tx_unit.get(
                            "cds_ranges",
                            [],
                        )

                        if not cds_ranges:
                            skipped_transcript_no_cds += 1
                            print(
                                f"[SKIP] {gene_id} {tx_id}: "
                                "transcript exists but has no CDS"
                            )
                            continue

                        total_transcripts_with_cds += 1

                        cds_chr = (
                            tx_unit.get("chr")
                            or chrom
                        )
                        cds_strand = (
                            tx_unit.get("strand")
                            or strand
                        )

                        cds_start = min(
                            s for s, e in cds_ranges
                        )
                        cds_end = max(
                            e for s, e in cds_ranges
                        )

                        cds_ranges_array = (
                            make_cds_ranges_array(
                                cds_ranges,
                                cds_strand,
                            )
                        )

                        region = extract_cds_variants(
                            vcf,
                            cds_chr,
                            cds_ranges,
                        )

                        dels_chr = dels_by_chr.get(
                            cds_chr,
                            [],
                        )

                        extra_dels = (
                            extract_cds_overlapping_dels_from_list(
                                dels_chr,
                                cds_ranges,
                            )
                        )

                        var_map = {
                            (
                                s["pos"],
                                s["ref"],
                                s["alt"],
                            ): s
                            for s in region
                        }

                        for s in extra_dels:
                            k = (
                                s["pos"],
                                s["ref"],
                                s["alt"],
                            )

                            if k in var_map:
                                continue

                            geno = (
                                fetch_site_genotypes_from_vcf(
                                    vcf,
                                    cds_chr,
                                    s["pos"],
                                    s["ref"],
                                    s["alt"],
                                )
                            )

                            if geno is None:
                                continue

                            s2 = dict(s)
                            s2["genotypes"] = geno
                            var_map[k] = s2

                        region = [
                            var_map[k]
                            for k in sorted(
                                var_map,
                                key=lambda x: x[0],
                            )
                        ]

                        if not region:
                            print(
                                f"[SKIP] {gene_id} {tx_id}: "
                                "no variants found in CDS "
                                f"{cds_chr}:{cds_start}-{cds_end}"
                            )
                            continue

                        block_chrom = cds_chr
                        block_start_1based = cds_start
                        block_end_1based = cds_end
                        block_id = (
                            f"{block_chrom}:"
                            f"{block_start_1based}-"
                            f"{block_end_1based}|"
                            f"{gene_id}|"
                            f"{tx_id}"
                        )

                        finalize_and_write_block(
                            out=out,
                            block_chrom=block_chrom,
                            block_start_1based=(
                                block_start_1based
                            ),
                            block_end_1based=(
                                block_end_1based
                            ),
                            block_id=block_id,
                            gene_id=gene_id,
                            strand=cds_strand,
                            block_mode="CDS",
                            flank_value=".",
                            cds_ranges_array=(
                                cds_ranges_array
                            ),
                            region=region,
                            transcript_id=tx_id,
                        )
    finally:
        vcf.close()

    print(f"[INFO] Written to {output_path}")

    with open(
        summary_path,
        "w",
        encoding="utf-8",
    ) as s:
        s.write("metric\tvalue\n")
        s.write(f"gene_count\t{gene_count}\n")
        s.write(
            f"total_genes_in_gff\t"
            f"{total_genes_in_gff}\n"
        )
        s.write(
            f"gene_id_list_raw_count\t"
            f"{len(raw_gene_ids) if raw_gene_ids is not None else 0}\n"
        )
        s.write("annotation_format\tGFF3\n")
        s.write(
            "transcript_selection_mode\t"
            "all_transcripts\n"
        )
        s.write(
            f"min_sample_count\t"
            f"{config.min_sample_count}\n"
        )
        s.write(f"mode\t{config.mode}\n")
        s.write("format_version\t3\n")
        s.write(
            "coordinate_system\t"
            "1-based-closed\n"
        )
        s.write(
            f"meta_query_region\t"
            f"{META_QUERY_REGION}\n"
        )
        s.write(
            f"total_transcripts_seen\t"
            f"{total_transcripts_seen}\n"
        )
        s.write(
            f"total_transcripts_with_cds\t"
            f"{total_transcripts_with_cds}\n"
        )
        s.write(
            f"skipped_transcript_no_cds\t"
            f"{skipped_transcript_no_cds}\n"
        )
        s.write(
            f"total_blocks_written\t"
            f"{total_blocks_written}\n"
        )
        s.write(
            f"total_transcript_blocks_written\t"
            f"{total_transcript_blocks_written}\n"
        )
        s.write(
            f"total_haplotypes_filtered\t"
            f"{total_hap_count}\n"
        )

        if total_blocks_written > 0:
            s.write(
                "avg_haplotypes_per_block_filtered\t"
                f"{total_hap_count / total_blocks_written:.2f}\n"
            )
        else:
            s.write(
                "avg_haplotypes_per_block_filtered\t0\n"
            )

        s.write(
            "\n# sample_count\thaplotype_count\n"
        )

        for sc in sorted(
            hap_sample_count_map.keys()
        ):
            s.write(
                f"{sc}\t"
                f"{hap_sample_count_map[sc]}\n"
            )

    print(
        f"[INFO] Summary written to {summary_path}"
    )

    # Large haplotype TSV files are automatically sorted, bgzip-compressed,
    # and Tabix-indexed. The original TSV is kept unchanged.
    compress_large_haplotype_tsv_if_needed(output_path)

    return {
        "output_path": output_path,
        "summary_path": summary_path,
        "all_del_path": del_position_file,
        "gene_count": gene_count,
        "total_blocks_written": total_blocks_written,
    }