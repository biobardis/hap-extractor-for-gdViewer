"""Generic GFF3 validation, parsing, and transcript normalization."""

import re
from collections import defaultdict, deque
from urllib.parse import unquote

from .constants import TRANSCRIPT_TYPES
from .intervals import merge_intervals

def clean_feature_type(ftype):
    return str(ftype).strip().lower()

def parse_gff_attributes(attr_str):
    """
    Parse standard GFF3 column 9:
        key=value;key=value

    Values are URL-decoded.
    Unknown attributes are preserved.
    """
    attrs = {}

    for item in attr_str.strip().split(";"):
        item = item.strip()
        if not item:
            continue

        if "=" not in item:
            # Non-standard GFF3 fragment. Keep it rather than crashing.
            attrs[item] = ""
            continue

        key, value = item.split("=", 1)
        key = key.strip()
        value = unquote(value.strip())
        attrs[key] = value

    return attrs

def split_parent_values(parent_value):
    if parent_value is None:
        return []

    return [
        x.strip()
        for x in str(parent_value).split(",")
        if x.strip()
    ]

def validate_gff3(gff_path, max_data_lines=5000):
    """
    Validate that the input appears to be GFF3 rather than GTF.

    Rules:
    - If ##gff-version 3 exists -> strong confirmation.
    - Otherwise inspect data lines:
      * exactly >= 9 tab-separated columns
      * GFF3-style key=value attributes
      * reject obvious GTF-style gene_id "..." / transcript_id "..."
    """
    found_header = False
    data_lines = 0
    gff3_attr_evidence = 0
    gtf_attr_evidence = 0
    malformed = 0

    with open(gff_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("##gff-version"):
                if line.strip().lower() == "##gff-version 3":
                    found_header = True

            if line.startswith("#") or not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                malformed += 1
                continue

            data_lines += 1
            attrs = cols[8]

            if re.search(r'\b(?:gene_id|transcript_id)\s+"[^"]+"', attrs):
                gtf_attr_evidence += 1

            if re.search(r'(?:^|;)\s*[A-Za-z0-9_.:-]+\s*=', attrs):
                gff3_attr_evidence += 1

            if data_lines >= max_data_lines:
                break

    if data_lines == 0:
        raise ValueError("No valid data lines were found in the annotation file.")

    if gtf_attr_evidence > 0 and gff3_attr_evidence == 0:
        raise ValueError(
            "Input appears to be GTF, but this script only supports GFF3."
        )

    if gff3_attr_evidence == 0:
        raise ValueError(
            "Column 9 does not look like GFF3 key=value attributes. "
            "This script only supports GFF3."
        )

    if found_header:
        print("[GFF3] Confirmed by '##gff-version 3' header.")
    else:
        print(
            "[WARNING] Missing '##gff-version 3' header, "
            "but the file structure appears to be GFF3. Continuing."
        )

    if malformed:
        print(
            f"[WARNING] Found {malformed} malformed records with fewer than 9 columns "
            f"among the inspected portion."
        )

def parse_all_gff3_features(gff_path):
    """
    First pass: read all GFF3 records without assuming parent-before-child order.

    Returns:
        features: list[dict]
        feature_by_id: id -> first feature carrying that ID
        children_by_parent: parent_id -> list[feature]
        duplicate_id_count
    """
    features = []
    feature_by_id = {}
    children_by_parent = defaultdict(list)
    duplicate_id_count = 0

    order = 0

    with open(gff_path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            if line.startswith("#") or not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            chrom = cols[0]
            source = cols[1]
            ftype = cols[2]

            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                print(f"[WARNING] Skip line {line_no}: invalid start/end.")
                continue

            if end < start:
                start, end = end, start

            score = cols[5]
            strand = cols[6]
            phase = cols[7]
            attrs = parse_gff_attributes(cols[8])

            feature_id = attrs.get("ID", "").strip()
            parents = split_parent_values(attrs.get("Parent"))

            feat = {
                "order": order,
                "line_no": line_no,
                "chr": chrom,
                "source": source,
                "type": ftype,
                "type_norm": clean_feature_type(ftype),
                "start": start,
                "end": end,
                "score": score,
                "strand": strand,
                "phase": phase,
                "id": feature_id,
                "parents": parents,
                "attributes": attrs,
            }

            features.append(feat)
            order += 1

            if feature_id:
                if feature_id not in feature_by_id:
                    feature_by_id[feature_id] = feat
                else:
                    duplicate_id_count += 1

            for parent_id in parents:
                children_by_parent[parent_id].append(feat)

    return features, feature_by_id, children_by_parent, duplicate_id_count

def direct_has_exon_or_cds_child(feature_id, children_by_parent):
    for child in children_by_parent.get(feature_id, []):
        if child["type_norm"] in {"exon", "cds"}:
            return True
    return False

def is_explicit_gene_feature(feat):
    """
    Gene-level detection.

    - `gene` is treated as gene-level.
    - A top-level `pseudogene` (no Parent) is treated as gene-level.
    - A `pseudogene` carrying Parent is NOT forced to gene-level, because some
      Ensembl-style GFF3 files use `pseudogene` as a transcript-level feature.
    """
    if feat["type_norm"] == "gene":
        return True

    if feat["type_norm"] == "pseudogene" and not feat.get("parents"):
        return True

    return False

def is_transcript_candidate(feat, children_by_parent):
    """
    Transcript detection:
    1. known transcript-like type
    OR
    2. unknown/non-gene feature with ID that directly owns exon/CDS child

    Explicit gene/pseudogene is never classified as a transcript.
    """
    if not feat["id"]:
        return False

    if is_explicit_gene_feature(feat):
        return False

    if feat["type_norm"] in TRANSCRIPT_TYPES:
        return True

    # Some Ensembl-style annotations can use `pseudogene` for a
    # transcript-level feature under a gene.
    if feat["type_norm"] == "pseudogene" and feat.get("parents"):
        return True

    return direct_has_exon_or_cds_child(feat["id"], children_by_parent)

def find_ancestor_gene_ids(
    start_feature,
    feature_by_id,
    explicit_gene_ids,
    inferred_gene_parent_ids,
    max_depth=50,
):
    """
    Follow Parent links upward and return nearest gene-level ancestor IDs.

    Gene-level nodes:
    - explicit gene/pseudogene features
    - IDs inferred as gene parents of transcript candidates
    """
    found = []
    seen = set()
    q = deque()

    for p in start_feature.get("parents", []):
        q.append((p, 1))

    best_depth = None

    while q:
        pid, depth = q.popleft()

        if pid in seen or depth > max_depth:
            continue
        seen.add(pid)

        if best_depth is not None and depth > best_depth:
            continue

        if pid in explicit_gene_ids or pid in inferred_gene_parent_ids:
            found.append(pid)
            best_depth = depth
            continue

        parent_feat = feature_by_id.get(pid)
        if parent_feat is None:
            continue

        for pp in parent_feat.get("parents", []):
            q.append((pp, depth + 1))

    return list(dict.fromkeys(found))

def collect_descendant_ranges(
    root_id,
    children_by_parent,
    wanted_type,
    max_nodes=100000,
):
    """
    Recursively collect descendant exon/CDS ranges under a transcript.
    Recursive traversal makes the parser tolerant of an extra intermediate level.

    Because gene nodes are not passed here, this does not cause a gene to become
    a transcript merely because it has deep CDS descendants.
    """
    out = []
    q = deque(children_by_parent.get(root_id, []))
    visited_keys = set()
    n = 0

    while q:
        feat = q.popleft()
        n += 1
        if n > max_nodes:
            raise RuntimeError(
                f"Too many descendants while traversing {root_id}; possible malformed cycle."
            )

        key = (
            feat["line_no"],
            feat["id"],
            feat["type_norm"],
            feat["start"],
            feat["end"],
        )
        if key in visited_keys:
            continue
        visited_keys.add(key)

        if feat["type_norm"] == wanted_type:
            out.append((feat["start"], feat["end"]))

        if feat["id"]:
            for child in children_by_parent.get(feat["id"], []):
                q.append(child)

    return merge_intervals(out)

def make_transcript_unit(
    gene_id,
    tx_feat,
    children_by_parent,
    synthetic=False,
    synthetic_cds_ranges=None,
    synthetic_exon_ranges=None,
):
    if synthetic:
        cds_ranges = merge_intervals(synthetic_cds_ranges or [])
        exon_ranges = merge_intervals(synthetic_exon_ranges or [])

        return {
            "gene_id": gene_id,
            "transcript_id": tx_feat["id"],
            "chr": tx_feat["chr"],
            "strand": tx_feat["strand"],
            "tx_start": tx_feat["start"],
            "tx_end": tx_feat["end"],
            "cds_ranges": cds_ranges,
            "exon_ranges": exon_ranges,
            "order": tx_feat["order"],
        }

    cds_ranges = collect_descendant_ranges(
        tx_feat["id"], children_by_parent, "cds"
    )
    exon_ranges = collect_descendant_ranges(
        tx_feat["id"], children_by_parent, "exon"
    )

    return {
        "gene_id": gene_id,
        "transcript_id": tx_feat["id"],
        "chr": tx_feat["chr"],
        "strand": tx_feat["strand"],
        "tx_start": tx_feat["start"],
        "tx_end": tx_feat["end"],
        "cds_ranges": cds_ranges,
        "exon_ranges": exon_ranges,
        "order": tx_feat["order"],
    }

def load_gene_regions(gff_path):
    """
    Generic GFF3 normalization layer.

    Returns:
        genes:
            [(gene_id, chrom, start, end, strand), ...]

        transcript_units_map:
            gene_id -> [all transcript units, including no-CDS transcripts]

        exon_units_map:
            same transcript units filtered to those having exon ranges

    Important behavior:
    - No ID-prefix assumptions.
    - No protein_coding/basic/canonical filtering.
    - All transcript candidates are retained.
    - Direct gene -> exon/CDS children create a synthetic transcript.
    """
    validate_gff3(gff_path)

    (
        features,
        feature_by_id,
        children_by_parent,
        duplicate_id_count,
    ) = parse_all_gff3_features(gff_path)

    if duplicate_id_count:
        print(
            f"[WARNING] {duplicate_id_count} duplicate feature IDs were encountered. "
            "The first occurrence is used as the parent lookup target."
        )

    explicit_gene_ids = {
        feat["id"]
        for feat in features
        if feat["id"] and is_explicit_gene_feature(feat)
    }

    transcript_candidates = [
        feat
        for feat in features
        if is_transcript_candidate(feat, children_by_parent)
    ]
    transcript_candidate_ids = {x["id"] for x in transcript_candidates}

    # Any direct Parent of a transcript candidate can act as a gene-level parent.
    # This allows unusual but still structurally meaningful GFF3 styles.
    inferred_gene_parent_ids = set()
    for tx in transcript_candidates:
        for pid in tx["parents"]:
            if pid not in transcript_candidate_ids:
                inferred_gene_parent_ids.add(pid)

    gene_ids = set(explicit_gene_ids) | inferred_gene_parent_ids

    # Transcript -> gene assignment.
    tx_gene_map = defaultdict(list)
    unassigned_transcripts = []

    for tx in transcript_candidates:
        gene_ancestors = find_ancestor_gene_ids(
            tx,
            feature_by_id,
            explicit_gene_ids,
            inferred_gene_parent_ids,
        )

        if not gene_ancestors:
            # Last-resort structural fallback:
            # if a transcript has a Parent, treat the first non-transcript parent as gene.
            fallback = [
                p for p in tx["parents"]
                if p not in transcript_candidate_ids
            ]
            if fallback:
                gene_ancestors = [fallback[0]]
                gene_ids.add(fallback[0])

        if not gene_ancestors:
            unassigned_transcripts.append(tx["id"])
            continue

        for gid in gene_ancestors:
            tx_gene_map[gid].append(tx)

    # Build transcript units.
    transcript_units_map = defaultdict(list)

    for gene_id, tx_list in tx_gene_map.items():
        seen_tx = set()

        for tx in sorted(tx_list, key=lambda x: x["order"]):
            if tx["id"] in seen_tx:
                continue
            seen_tx.add(tx["id"])

            unit = make_transcript_unit(
                gene_id=gene_id,
                tx_feat=tx,
                children_by_parent=children_by_parent,
                synthetic=False,
            )
            transcript_units_map[gene_id].append(unit)

    # Direct gene -> exon/CDS children become one synthetic transcript per gene.
    synthetic_count = 0

    def gene_order_key(gene_id):
        gene_feat = feature_by_id.get(gene_id)
        if gene_feat is not None:
            return (gene_feat["order"], gene_id)

        direct_children = children_by_parent.get(gene_id, [])
        if direct_children:
            return (min(x["order"] for x in direct_children), gene_id)

        return (10**18, gene_id)

    for gene_id in sorted(gene_ids, key=gene_order_key):
        direct_children = children_by_parent.get(gene_id, [])

        direct_exons = [
            (x["start"], x["end"])
            for x in direct_children
            if x["type_norm"] == "exon"
        ]
        direct_cds = [
            (x["start"], x["end"])
            for x in direct_children
            if x["type_norm"] == "cds"
        ]

        if not direct_exons and not direct_cds:
            continue

        gene_feat = feature_by_id.get(gene_id)

        child_spans = direct_exons + direct_cds
        child_start = min(s for s, e in child_spans)
        child_end = max(e for s, e in child_spans)

        if gene_feat is not None:
            chrom = gene_feat["chr"]
            strand = gene_feat["strand"]
            tx_start = gene_feat["start"]
            tx_end = gene_feat["end"]
            order = gene_feat["order"]
            attrs = gene_feat["attributes"]
        else:
            first_child = direct_children[0]
            chrom = first_child["chr"]
            strand = first_child["strand"]
            tx_start = child_start
            tx_end = child_end
            order = min(x["order"] for x in direct_children)
            attrs = {}

        synthetic_feat = {
            "id": f"transcript{synthetic_count + 1}",
            "chr": chrom,
            "strand": strand,
            "start": tx_start,
            "end": tx_end,
            "order": order,
            "attributes": attrs,
        }

        unit = make_transcript_unit(
            gene_id=gene_id,
            tx_feat=synthetic_feat,
            children_by_parent=children_by_parent,
            synthetic=True,
            synthetic_cds_ranges=direct_cds,
            synthetic_exon_ranges=direct_exons,
        )

        transcript_units_map[gene_id].append(unit)
        synthetic_count += 1

    # Sort all transcripts by original GFF order.
    for gid in transcript_units_map:
        transcript_units_map[gid].sort(
            key=lambda x: (x.get("order", 10**18), x["transcript_id"])
        )

    # Build normalized gene coordinates.
    normalized_genes = []
    inferred_gene_count = 0

    for gene_id in gene_ids:
        gene_feat = feature_by_id.get(gene_id)

        if gene_feat is not None:
            chrom = gene_feat["chr"]
            start = gene_feat["start"]
            end = gene_feat["end"]
            strand = gene_feat["strand"]
            order = gene_feat["order"]
        else:
            units = transcript_units_map.get(gene_id, [])
            if not units:
                continue

            chroms = [u["chr"] for u in units if u.get("chr")]
            if not chroms:
                continue
            chrom = chroms[0]

            starts = [u["tx_start"] for u in units if u.get("tx_start") is not None]
            ends = [u["tx_end"] for u in units if u.get("tx_end") is not None]
            if not starts or not ends:
                continue

            start = min(starts)
            end = max(ends)
            strands = [u["strand"] for u in units if u.get("strand")]
            strand = strands[0] if strands else "."
            order = min(u.get("order", 10**18) for u in units)
            inferred_gene_count += 1

        normalized_genes.append(
            (order, gene_id, chrom, start, end, strand)
        )

    normalized_genes.sort(key=lambda x: (x[0], x[2], x[3], x[1]))
    genes = [
        (gene_id, chrom, start, end, strand)
        for order, gene_id, chrom, start, end, strand in normalized_genes
    ]

    exon_units_map = defaultdict(list)
    for gid, units in transcript_units_map.items():
        exon_units_map[gid] = [
            u for u in units if u.get("exon_ranges")
        ]

    print("[GFF3] Parsed features:", len(features))
    print("[GFF3] Explicit gene/pseudogene IDs:", len(explicit_gene_ids))
    print("[GFF3] Normalized gene IDs:", len(genes))
    print("[GFF3] Transcript candidates:", len(transcript_candidates))
    print("[GFF3] Generated transcripts:", synthetic_count)
    print("[GFF3] Inferred gene records:", inferred_gene_count)
    print(
        "[GFF3] Total normalized transcript units:",
        sum(len(v) for v in transcript_units_map.values()),
    )

    if unassigned_transcripts:
        print(
            f"[WARNING] {len(unassigned_transcripts)} transcript candidates "
            "could not be assigned to any gene and were skipped."
        )
        for txid in unassigned_transcripts[:10]:
            print("  [UNASSIGNED]", txid)
        if len(unassigned_transcripts) > 10:
            print("  ...")

    return genes, transcript_units_map, exon_units_map