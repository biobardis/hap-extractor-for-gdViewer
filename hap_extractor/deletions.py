"""Deletion extraction and spanning-deletion recovery."""

import os
import sys

import pysam

from .intervals import any_overlap, intervals_overlap, point_in_intervals
from .paths import ensure_parent_dir

def is_biallelic_del(rec):
    """
    Return True only for a single-ALT deletion record.

    Kept:
    1. sequence-resolved deletion: len(REF) > len(ALT)
    2. symbolic deletion: ALT=<DEL>

    Excluded:
    - records without ALT
    - multiallelic records
    - *, <*>, .
    """
    if not rec.alts or len(rec.alts) != 1:
        return False

    alt = rec.alts[0]
    if alt in ("<*>", "*", "."):
        return False

    if alt == "<DEL>":
        return True

    return len(rec.ref) > len(alt)

def is_multiallelic_del(rec):
    """Statistics only: whether a multiallelic record contains a deletion allele."""
    if not rec.alts or len(rec.alts) <= 1:
        return False

    for alt in rec.alts:
        if alt in ("<*>", "*", "."):
            continue
        if alt == "<DEL>":
            return True
        if len(rec.ref) > len(alt):
            return True

    return False

def get_del_end_1based(rec):
    """Return a 1-based inclusive deletion end coordinate."""
    alt = rec.alts[0] if rec.alts else None

    if alt == "<DEL>":
        try:
            end = rec.info.get("END")
            if end is not None:
                return int(end)
        except Exception:
            pass

        # pysam rec.stop is also 1-based inclusive from the VCF perspective.
        try:
            if rec.stop is not None:
                return int(rec.stop)
        except Exception:
            pass

    return rec.pos + len(rec.ref) - 1

def get_del(vcf_path, output_path, report_every=100000):
    """
    Generate the all-DEL TSV used by the spanning-deletion recovery step.

    Only biallelic/single-ALT deletions are written. Multiallelic records that
    contain a deletion are counted but not written.
    """
    ensure_parent_dir(output_path)

    print(f"[GET_DEL] No --all-del supplied; generating: {output_path}")
    print(f"[GET_DEL] Source VCF: {vcf_path}")

    vcf = pysam.VariantFile(vcf_path)
    total = 0
    del_count = 0
    skipped_multiallelic_del_count = 0

    try:
        with open(output_path, "w", encoding="utf-8") as out:
            out.write("chrom\tstart\tend\tref\talt\n")

            # Sequential iteration is sufficient for generating the whole-file list.
            for rec in vcf:
                total += 1

                if total % report_every == 0:
                    print(
                        f"[GET_DEL] processed {total:,} records, "
                        f"written {del_count:,} biallelic DELs, "
                        f"skipped {skipped_multiallelic_del_count:,} "
                        "multiallelic DEL records",
                        file=sys.stderr,
                    )

                if is_multiallelic_del(rec):
                    skipped_multiallelic_del_count += 1
                    continue

                if not is_biallelic_del(rec):
                    continue

                alt = rec.alts[0]
                end = get_del_end_1based(rec)

                out.write(
                    f"{rec.chrom}\t{rec.pos}\t{end}\t{rec.ref}\t{alt}\n"
                )
                del_count += 1
    finally:
        vcf.close()

    print(
        f"[GET_DEL DONE] total records: {total:,}; "
        f"biallelic DELs written: {del_count:,}; "
        f"multiallelic DEL records skipped: {skipped_multiallelic_del_count:,}"
    )

    return output_path

def load_all_dels_tsv(del_tsv_path):
    """
    Read:
        chrom start end ref alt

    Coordinates: 1-based inclusive.

    If path is None / empty / missing, return an empty dict.
    """
    if not del_tsv_path:
        print("[INFO] DEL_POSITION_FILE disabled.")
        return {}

    if not os.path.exists(del_tsv_path):
        print(
            f"[WARNING] DEL_POSITION_FILE not found: {del_tsv_path}. "
            "Extra spanning-deletion recovery is disabled."
        )
        return {}

    dels_by_chr = {}

    with open(del_tsv_path, "r", encoding="utf-8") as f:
        header = f.readline()

        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                continue

            chrom, s, e, ref, alt = parts[:5]
            s = int(s)
            e = int(e)

            if e < s:
                s, e = e, s

            dels_by_chr.setdefault(chrom, []).append(
                (s, e, ref, alt)
            )

    for chrom in dels_by_chr:
        dels_by_chr[chrom].sort(key=lambda x: x[0])

    return dels_by_chr

def extract_cds_overlapping_dels_from_list(dels_chr, cds_ranges_1based):
    """
    Keep deletions whose START is outside CDS but whose interval overlaps CDS.
    """
    if not cds_ranges_1based or not dels_chr:
        return []

    cds_min = min(s for s, e in cds_ranges_1based)
    cds_max = max(e for s, e in cds_ranges_1based)

    out = []
    seen = set()

    for ds, de, ref, alt in dels_chr:
        if ds > cds_max:
            break
        if de < cds_min:
            continue
        if point_in_intervals(ds, cds_ranges_1based):
            continue
        if not any_overlap(ds, de, cds_ranges_1based):
            continue

        key = (ds, de, ref, alt)
        if key in seen:
            continue
        seen.add(key)

        out.append({
            "pos": ds,
            "end": de,
            "ref": ref,
            "alt": alt,
            "type": "DEL_CDSSPAN",
            "len": str(de - ds + 1),
        })

    return out

def extract_gene_overlapping_dels_from_list(
    dels_chr,
    region_start_1based,
    region_end_1based,
):
    """
    Keep deletions whose START is outside the gene block but overlap the block.
    """
    if not dels_chr:
        return []

    out = []
    seen = set()

    for ds, de, ref, alt in dels_chr:
        if ds > region_end_1based:
            break

        if de < region_start_1based:
            continue

        if region_start_1based <= ds <= region_end_1based:
            continue

        if not intervals_overlap(
            ds,
            de,
            region_start_1based,
            region_end_1based,
        ):
            continue

        key = (ds, de, ref, alt)
        if key in seen:
            continue
        seen.add(key)

        out.append({
            "pos": ds,
            "end": de,
            "ref": ref,
            "alt": alt,
            "type": "DEL_GENESPAN",
            "len": str(de - ds + 1),
        })

    return out