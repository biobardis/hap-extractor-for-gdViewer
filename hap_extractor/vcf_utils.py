"""VCF genotype and variant extraction helpers."""

from .intervals import any_overlap, point_in_intervals

def get_diploid_gt_strings(call):
    """
    Return exactly two allele strings.

    Missing/haploid GTs are padded with ".".
    Polyploid GTs are truncated to the first two alleles because the downstream
    haplotype model is diploid.
    """
    gt = call.get("GT")
    if gt is None:
        return None

    alleles = [
        "." if a is None else str(a)
        for a in gt
    ]

    if len(alleles) == 0:
        return None
    if len(alleles) == 1:
        alleles.append(".")
    if len(alleles) > 2:
        alleles = alleles[:2]

    return tuple(alleles)

def fetch_site_genotypes_from_vcf(vcf, chrom, pos_1based, ref, alt):
    for rec in vcf.fetch(chrom, pos_1based - 1, pos_1based):
        if rec.alts is None:
            continue

        rec_alt = rec.alts[0]
        if (
            rec.pos == pos_1based
            and rec.ref == ref
            and rec_alt == alt
        ):
            geno = {}
            for s in vcf.header.samples:
                gt = get_diploid_gt_strings(rec.samples[s])
                if gt is None:
                    continue
                geno[s] = gt
            return geno

    return None

def get_variant_kind(ref, alt):
    if alt is None or alt == ".":
        return "UNKNOWN"

    if alt == "<DEL>":
        return "DEL"

    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(ref) < len(alt):
        return "INS"
    if len(ref) > len(alt):
        return "DEL"

    return "MNP_OR_COMPLEX"

def insertion_anchor_inside_intervals(pos_1based, intervals_1based):
    """
    VCF insertion POS is the left anchor base.

    For CDS [s,e]:
      insertion after e is outside;
      insertion after s-1 is outside;
      only s <= anchor < e is inside that interval.
    """
    for s, e in intervals_1based:
        if s <= pos_1based < e:
            return True
    return False

def variant_overlaps_intervals_strict(
    pos_1based,
    ref,
    alt,
    intervals_1based,
):
    """
    Determine whether the changed portion of REF/ALT truly affects target intervals.
    """
    if alt is None or alt == ".":
        return point_in_intervals(pos_1based, intervals_1based)

    prefix = 0
    max_prefix = min(len(ref), len(alt))

    while prefix < max_prefix and ref[prefix] == alt[prefix]:
        prefix += 1

    ref_suffix = len(ref)
    alt_suffix = len(alt)

    while (
        ref_suffix > prefix
        and alt_suffix > prefix
        and ref[ref_suffix - 1] == alt[alt_suffix - 1]
    ):
        ref_suffix -= 1
        alt_suffix -= 1

    ref_changed_len = ref_suffix - prefix
    alt_changed_len = alt_suffix - prefix

    if ref_changed_len == 0 and alt_changed_len > 0:
        insertion_anchor = pos_1based + prefix - 1
        return insertion_anchor_inside_intervals(
            insertion_anchor,
            intervals_1based,
        )

    if ref_changed_len > 0:
        affected_start = pos_1based + prefix
        affected_end = pos_1based + ref_suffix - 1
        return any_overlap(
            affected_start,
            affected_end,
            intervals_1based,
        )

    return False

def extract_region_variants(vcf, chrom, start0, end0):
    """
    Fetch a continuous 0-based half-open interval from pysam.
    Returned positions remain VCF 1-based positions.
    """
    region = []

    for rec in vcf.fetch(chrom, start0, end0):
        # Preserve the original GENE-mode behavior: use the first ALT.
        if rec.alts is None:
            continue

        pos = rec.pos
        ref = rec.ref
        alt = rec.alts[0]

        vtype = get_variant_kind(ref, alt)

        if alt == "<DEL>":
            try:
                end = rec.info.get("END")
            except Exception:
                end = None
            if end is None:
                try:
                    end = rec.stop
                except Exception:
                    end = None
            if end is not None:
                vlen = str(max(1, int(end) - pos + 1))
            else:
                vlen = str(len(ref))
        elif vtype == "SNP":
            vlen = "."
        elif vtype == "INS":
            vlen = str(len(alt))
        else:
            vlen = str(len(ref))

        geno = {}
        for s in vcf.header.samples:
            gt = get_diploid_gt_strings(rec.samples[s])
            if gt is None:
                continue
            geno[s] = gt

        region.append({
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "type": vtype,
            "len": vlen,
            "genotypes": geno,
        })

    return region

def extract_cds_variants(vcf, chrom, cds_ranges_1based):
    """
    Fetch every CDS interval, de-duplicate records, and retain only variants whose
    changed sequence truly affects the CDS union.
    """
    if not cds_ranges_1based:
        return []

    var_map = {}

    for s1, e1 in cds_ranges_1based:
        start0 = max(0, s1 - 1)
        end0 = e1

        for rec in vcf.fetch(chrom, start0, end0):
            if rec.alts is None or len(rec.alts) != 1:
                continue

            pos = rec.pos
            ref = rec.ref
            alt = rec.alts[0]

            if not variant_overlaps_intervals_strict(
                pos,
                ref,
                alt,
                cds_ranges_1based,
            ):
                continue

            key = (pos, ref, alt)
            if key in var_map:
                continue

            vtype = get_variant_kind(ref, alt)

            if alt == "<DEL>":
                try:
                    end = rec.info.get("END")
                except Exception:
                    end = None
                if end is None:
                    try:
                        end = rec.stop
                    except Exception:
                        end = None
                if end is not None:
                    vlen = str(max(1, int(end) - pos + 1))
                else:
                    vlen = str(len(ref))
            elif vtype == "SNP":
                vlen = "."
            elif vtype == "INS":
                vlen = str(len(alt))
            else:
                vlen = str(len(ref))

            geno = {}
            for s in vcf.header.samples:
                gt = get_diploid_gt_strings(rec.samples[s])
                if gt is None:
                    continue
                geno[s] = gt

            var_map[key] = {
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "type": vtype,
                "len": vlen,
                "genotypes": geno,
            }

    return [
        var_map[k]
        for k in sorted(var_map, key=lambda x: x[0])
    ]