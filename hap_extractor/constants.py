"""Shared constants for the haplotype extractor."""

META_CHROM = "META"
META_POS = 1
META_QUERY_REGION = "META:1-1"

TRANSCRIPT_TYPES = {
    "transcript",
    "mrna",
    "pseudogenic_transcript",
    "ncrna",
    "lncrna",
    "lnc_rna",
    "rrna",
    "trna",
    "mirna",
    "snrna",
    "snorna",
    "scrna",
    "srp_rna",
    "rnase_mrp_rna",
    "rnase_p_rna",
    "antisense_rna",
    "primary_transcript",
    "nmd_transcript_variant",
    "protein_coding_transcript",
}

HAPBLOCK_TSV_HEADER = [
    "chr",
    "start",
    "end",
    "record_type",
    "block_id",
    "anchor_id",
    "value_json",
]