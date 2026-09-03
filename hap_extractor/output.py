"""Output-format helpers."""

import json

from .constants import HAPBLOCK_TSV_HEADER, META_CHROM, META_POS, META_QUERY_REGION

def compact_json(obj):
    return json.dumps(
        obj,
        ensure_ascii=False,
        separators=(",", ":"),
    )

def clean_tsv_value(x):
    if x is None:
        return "."

    s = str(x)
    s = s.replace("\t", " ")
    s = s.replace("\n", " ")
    s = s.replace("\r", " ")
    return s

def write_tsv_row(out, row):
    out.write(
        "\t".join(clean_tsv_value(x) for x in row)
        + "\n"
    )

def make_variant_array(region, key):
    return [s[key] for s in region]

def make_cds_ranges_array(cds_ranges, strand):
    """
    Positive strand: genomic order.
    Negative strand: reverse to transcript direction.
    """
    ranges = (
        cds_ranges
        if strand == "+"
        else cds_ranges[::-1]
    )
    return [[s, e] for s, e in ranges]

def write_format_header(out):
    out.write("#coordinate_system: 1-based-closed\n")
    out.write(f"#meta_query_region: {META_QUERY_REGION}\n")
    out.write("#" + "\t".join(HAPBLOCK_TSV_HEADER) + "\n")


def write_meta_rows(out, total_samples, samples, min_sample_count, mode, gene_id_list):
    meta_payload = {
        "total_samples": total_samples,
        "min_sample_count": min_sample_count,
        "mode": mode,
        "gene_id_list_count": len(gene_id_list) if gene_id_list else 0,
        "annotation_format": "GFF3",
        "transcript_selection_mode": "all_transcripts",
        "samples": samples,
    }

    write_tsv_row(out, [
        META_CHROM,
        META_POS,
        META_POS,
        "META",
        "GLOBAL",
        ".",
        compact_json(meta_payload),
    ])

def write_block_and_hap_rows(
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
    hap_rows,
):
    block_payload = {
        "strand": strand,
        "mode": block_mode,
        "flank": flank_value,
        "cds_ranges": cds_ranges_array,
        "variant_positions": make_variant_array(region, "pos"),
        "variant_refs": make_variant_array(region, "ref"),
        "variant_alts": make_variant_array(region, "alt"),
        "variant_types": make_variant_array(region, "type"),
        "variant_lengths": make_variant_array(region, "len"),
    }

    write_tsv_row(out, [
        block_chrom,
        block_start_1based,
        block_end_1based,
        "BLOCK",
        block_id,
        gene_id,
        compact_json(block_payload),
    ])

    for hap in hap_rows:
        hap_payload = {
            "gene_id": gene_id,
            "hap_id": hap["hap_id"],
            "sample_count": hap["sample_count"],
            "samples": hap["samples"],
            "haplotype": hap["haplotype"],
        }

        write_tsv_row(out, [
            block_chrom,
            block_start_1based,
            block_end_1based,
            hap["record_type"],
            block_id,
            gene_id,
            compact_json(hap_payload),
        ])