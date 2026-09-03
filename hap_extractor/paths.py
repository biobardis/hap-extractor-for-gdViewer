"""Path and lightweight input helpers."""

import os


def load_gene_ids_from_list(gene_list_value):
    """Resolve --gene-list selection.

    ``all`` selects every gene and is represented internally by ``None``.
    Otherwise, ``gene_list_value`` must be a text file containing one gene ID
    per line. Empty/comment-only files are rejected by the CLI before this
    function is called.
    """
    if gene_list_value.lower() == "all":
        return None

    out = []

    with open(gene_list_value, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            out.append(line.split()[0])

    # Keep input order but remove duplicates.
    return list(dict.fromkeys(out))


def get_vcf_basename(vcf_path):
    """
    Get the input VCF basename without common VCF suffixes.

    Examples:
        /data/sample.vcf.gz  -> sample
        /data/sample.vcf.bgz -> sample
        /data/sample.vcf     -> sample
    """
    basename = os.path.basename(
        os.path.expanduser(vcf_path)
    )

    suffixes = (
        ".vcf.gz",
        ".vcf.bgz",
        ".vcf",
    )

    for suffix in suffixes:
        if basename.lower().endswith(suffix):
            return basename[:-len(suffix)]

    # Fallback for unusual filenames.
    return os.path.splitext(basename)[0]


def make_default_all_del_path(vcf_path, output_dir):
    """
    Generate the default all-DEL filename from the input VCF name.

    Example:
        input VCF:
            litchi_314acc_reseq_snp_indel_filter.vcf.gz

        output:
            all_dels.biallelic_only.litchi_314acc_reseq_snp_indel_filter.tsv
    """
    vcf_basename = get_vcf_basename(vcf_path)

    filename = (
        f"all_dels.biallelic_only.{vcf_basename}.tsv"
    )

    return os.path.join(
        os.path.expanduser(output_dir),
        filename,
    )


def make_default_output_names(mode, debug, min_sample_count, vcf_path):
    """Generate default main output and summary filenames."""
    mode_tag = mode.lower()
    debug_suffix = "_test" if debug else ""
    filtered_suffix = "_filtered" if min_sample_count > 1 else ""

    vcf_basename = get_vcf_basename(vcf_path)

    hap_output = (
        f"haps{filtered_suffix}.{mode_tag}.{debug_suffix}.{vcf_basename}.tsv"
    )

    summary_output = (
        f"summary{filtered_suffix}.{mode_tag}.{debug_suffix}.{vcf_basename}.tsv"
    )

    return hap_output, summary_output


def resolve_output_path(output_dir, value, default_name):
    """
    Resolve an output path.

    If value is not supplied:
        <output_dir>/<default_name>

    If value is a simple filename:
        <output_dir>/<value>

    If value contains a directory or is absolute:
        use it directly.
    """
    output_dir = os.path.expanduser(output_dir)

    if not value:
        return os.path.join(
            output_dir,
            default_name,
        )

    value = os.path.expanduser(value)

    if os.path.isabs(value) or os.path.dirname(value):
        return value

    return os.path.join(
        output_dir,
        value,
    )


def make_unique_path(path, related_suffixes=()):
    """
    Return a non-existing output path by appending ``_2``, ``_3``, ...
    before the file extension when needed.

    ``related_suffixes`` can be used for sidecar files that belong to the
    same output. For example, for a haplotype TSV, passing
    ``(".gz", ".gz.tbi")`` also treats ``example.tsv.gz`` and
    ``example.tsv.gz.tbi`` as collisions for ``example.tsv``.

    Examples:
        result.tsv      -> result.tsv
        result.tsv      -> result_2.tsv   (if result.tsv already exists)
        result_2.tsv    -> result_3.tsv   (if result.tsv and result_2.tsv exist)
    """
    path = os.path.expanduser(path)
    stem, ext = os.path.splitext(path)

    def is_occupied(candidate):
        if os.path.exists(candidate):
            return True

        return any(
            os.path.exists(candidate + suffix)
            for suffix in related_suffixes
        )

    if not is_occupied(path):
        return path

    index = 2
    while True:
        candidate = f"{stem}_{index}{ext}"
        if not is_occupied(candidate):
            return candidate
        index += 1


def ensure_parent_dir(path):
    """Create the parent directory if it does not exist."""
    parent = os.path.dirname(
        os.path.abspath(
            os.path.expanduser(path)
        )
    )

    if parent:
        os.makedirs(
            parent,
            exist_ok=True,
        )