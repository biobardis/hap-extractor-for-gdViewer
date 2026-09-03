"""Gene-list matching helpers."""

import re

def normalize_gene_id_for_match(gene_id):
    """
    Conservative convenience normalization for gene-list matching.

    - gene:ENSG... -> ENSG...
    - ENSG....17  -> ENSG....

    Does NOT strip arbitrary prefixes before ':'.
    """
    if gene_id is None:
        return ""

    x = str(gene_id).strip()
    if not x:
        return ""

    if x.startswith("gene:"):
        x = x[len("gene:"):]

    m = re.match(r"^(.*)\.(\d+)$", x)
    if m:
        x = m.group(1)

    return x

def build_gene_id_filter(gene_id_list):
    if gene_id_list is None:
        return None, None

    if isinstance(gene_id_list, str):
        items = gene_id_list.splitlines()
    else:
        items = list(gene_id_list)

    if not items:
        return None, None

    raw_gene_ids = set()
    match_gene_ids = set()

    for item in items:
        line = str(item).strip()
        if not line or line.startswith("#"):
            continue

        gid = line.split()[0]
        raw_gene_ids.add(gid)
        match_gene_ids.add(gid)

        norm = normalize_gene_id_for_match(gid)
        if norm:
            match_gene_ids.add(norm)

    if not raw_gene_ids:
        return None, None

    return raw_gene_ids, match_gene_ids

def gene_id_in_list(gene_id, match_gene_ids):
    if match_gene_ids is None:
        return True

    candidates = {str(gene_id)}
    norm = normalize_gene_id_for_match(gene_id)
    if norm:
        candidates.add(norm)

    return bool(candidates & match_gene_ids)