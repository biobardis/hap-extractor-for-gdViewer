"""Haplotype clustering and filtering."""

from collections import defaultdict

def cluster_all_haplotypes(region, samples):
    """
    Build two haplotypes per diploid sample from ordered GT allele indices.

    Example:
        site1 0|1
        site2 1|0

        hapA = 0,1
        hapB = 1,0
    """
    hap_dict = defaultdict(list)

    for sample in samples:
        hapA = []
        hapB = []

        for site in region:
            g = site["genotypes"]

            if sample not in g:
                hapA.append(".")
                hapB.append(".")
            else:
                a, b = g[sample]
                hapA.append(a)
                hapB.append(b)

        hapA_str = ",".join(hapA)
        hapB_str = ",".join(hapB)

        hap_dict[hapA_str].append(sample)
        hap_dict[hapB_str].append(sample)

    return hap_dict

def build_haplotype_rows(hap_dict, min_sample_count):
    hap_rows = []
    filtered_hap_samples = []
    kept_hap_count = 0
    hap_sample_count_updates = defaultdict(int)

    hid = 1

    for hap_str, sam_list in hap_dict.items():
        sc = len(sam_list)

        if sc < min_sample_count:
            filtered_hap_samples.extend(sam_list)
            continue

        hap_rows.append({
            "record_type": "HAP",
            "hap_id": f"HAP{hid}",
            "sample_count": sc,
            "samples": ",".join(sam_list),
            "haplotype": hap_str,
        })

        kept_hap_count += 1
        hap_sample_count_updates[sc] += 1
        hid += 1

    if filtered_hap_samples:
        hap_rows.append({
            "record_type": "HAP_FILTERED",
            "hap_id": "HAP_FILTERED",
            "sample_count": len(filtered_hap_samples),
            "samples": ",".join(filtered_hap_samples),
            "haplotype": ".",
        })

    return (
        hap_rows,
        kept_hap_count,
        hap_sample_count_updates,
    )