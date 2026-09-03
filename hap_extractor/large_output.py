"""Automatic sorting, bgzip compression, and Tabix indexing for large haplotype TSV files."""

import heapq
import os
import tempfile
import time

import pysam


LARGE_OUTPUT_THRESHOLD_BYTES = 100 * 1024 * 1024
LARGE_OUTPUT_CHUNK_LINES = 500000
LARGE_OUTPUT_PROGRESS_EVERY = 200000


# -------------------
# Helper functions
# -------------------
def human_num(n: int) -> str:
    return f"{n:,}"


def human_time(seconds: float) -> str:
    seconds = int(seconds)
    h = seconds // 3600
    m = (seconds % 3600) // 60
    s = seconds % 60

    if h > 0:
        return f"{h}h {m}m {s}s"
    if m > 0:
        return f"{m}m {s}s"
    return f"{s}s"


# -------------------
# External sorting
# -------------------
def sort_large_tsv(
    input_file,
    output_file,
    chunk_lines=LARGE_OUTPUT_CHUNK_LINES,
    progress_every=LARGE_OUTPUT_PROGRESS_EVERY,
    tmp_dir=None,
):
    """
    Sort a large TSV file by chrom and start.

    Temporary chunk files will be written to tmp_dir if provided.
    """

    temp_files = []
    header_lines = []
    total_lines = 0
    data_lines = 0
    chunk_idx = 0
    t0 = time.time()

    if tmp_dir is not None:
        os.makedirs(tmp_dir, exist_ok=True)

    def sort_key(line):
        cols = line.rstrip("\n").split("\t")
        chrom = cols[0]
        start = int(cols[1])
        return (chrom, start)

    def write_chunk(chunk, chunk_idx):
        chunk.sort(key=sort_key)

        tf = tempfile.NamedTemporaryFile(
            delete=False,
            mode="w",
            encoding="utf-8",
            prefix=f"chunk_{chunk_idx}_",
            suffix=".tsv",
            dir=tmp_dir,
        )

        tf.writelines(chunk)
        tf.close()

        return tf.name

    try:
        chunk = []

        with open(input_file, "r", encoding="utf-8") as f:
            for line in f:
                total_lines += 1

                if not line.strip():
                    continue

                if line.startswith("#"):
                    header_lines.append(line)
                    continue

                chunk.append(line)
                data_lines += 1

                if data_lines % progress_every == 0:
                    print(
                        f"[INFO] read {human_num(data_lines)} data lines; "
                        f"created {chunk_idx} temp chunks; "
                        f"elapsed {human_time(time.time() - t0)}"
                    )

                if len(chunk) >= chunk_lines:
                    chunk_idx += 1
                    tmp_path = write_chunk(chunk, chunk_idx)
                    temp_files.append(tmp_path)
                    print(f"[INFO] wrote temp chunk {chunk_idx}: {tmp_path}")
                    chunk = []

            if chunk:
                chunk_idx += 1
                tmp_path = write_chunk(chunk, chunk_idx)
                temp_files.append(tmp_path)
                print(f"[INFO] wrote temp chunk {chunk_idx}: {tmp_path}")

        print(
            f"[INFO] Finished splitting: "
            f"{human_num(data_lines)} data lines, "
            f"{chunk_idx} temp chunks"
        )

        # Merge sorted chunks.
        def file_iter(path):
            with open(path, "r", encoding="utf-8") as f:
                for line in f:
                    yield line

        def keyed_iter(it):
            for line in it:
                yield (sort_key(line), line)

        print(f"[INFO] Merging {len(temp_files)} temp chunks...")

        with open(output_file, "w", encoding="utf-8") as out:
            if header_lines:
                out.writelines(header_lines)

            iters = [file_iter(p) for p in temp_files]

            written = 0
            for _, line in heapq.merge(*(keyed_iter(it) for it in iters)):
                out.write(line)
                written += 1

                if written % progress_every == 0:
                    print(
                        f"[INFO] merged {human_num(written)} lines; "
                        f"elapsed {human_time(time.time() - t0)}"
                    )

        print(f"[INFO] Merge done: {output_file}")

    finally:
        # Always try to clean temporary chunk files.
        for p in temp_files:
            try:
                os.unlink(p)
            except FileNotFoundError:
                pass


def compress_large_haplotype_tsv_if_needed(
    input_file,
    threshold_bytes=LARGE_OUTPUT_THRESHOLD_BYTES,
    chunk_lines=LARGE_OUTPUT_CHUNK_LINES,
    progress_every=LARGE_OUTPUT_PROGRESS_EVERY,
):
    """
    Automatically create a sorted bgzip-compressed and Tabix-indexed copy
    when the haplotype TSV is larger than the configured threshold.

    The original TSV is kept unchanged.

    Outputs for ``example.tsv``:
        example.tsv.gz
        example.tsv.gz.tbi
    """

    input_file = os.path.abspath(os.path.expanduser(input_file))
    file_size = os.path.getsize(input_file)

    if file_size <= threshold_bytes:
        print(
            f"[INFO] Haplotype TSV size is "
            f"{file_size / (1024 * 1024):.2f} MB; "
            "automatic bgzip/Tabix compression is not required "
            "(threshold: >100 MB)."
        )
        return None

    gz_tsv = f"{input_file}.gz"
    tbi_path = f"{gz_tsv}.tbi"
    output_dir = os.path.dirname(input_file) or "."

    print(
        f"[INFO] Haplotype TSV size is "
        f"{file_size / (1024 * 1024):.2f} MB (>100 MB)."
    )
    print("[INFO] Starting automatic sort + bgzip + Tabix indexing...")

    start_t0 = time.time()

    # Keep temporary sorting data beside the output file. This avoids relying
    # on the system temporary partition for very large files.
    with tempfile.TemporaryDirectory(
        prefix="hap_tabix_sort_",
        dir=output_dir,
    ) as tmp_dir:
        sorted_tsv = os.path.join(
            tmp_dir,
            os.path.basename(input_file) + ".sorted.tsv",
        )

        print(f"[INFO] Temporary sorting directory: {tmp_dir}")
        print(f"[INFO] Sorting {input_file} -> {sorted_tsv} ...")

        sort_large_tsv(
            input_file,
            sorted_tsv,
            chunk_lines=chunk_lines,
            progress_every=progress_every,
            tmp_dir=tmp_dir,
        )

        print(
            f"[INFO] Sorting done in "
            f"{human_time(time.time() - start_t0)}"
        )

        print(f"[INFO] Compressing with bgzip -> {gz_tsv} ...")
        pysam.tabix_compress(
            sorted_tsv,
            gz_tsv,
            force=True,
        )
        print("[INFO] Done compression")

    print(f"[INFO] Building Tabix index -> {tbi_path} ...")

    pysam.tabix_index(
        gz_tsv,
        seq_col=0,
        start_col=1,
        end_col=2,
        meta_char="#",
        force=True,
        preset=None,
    )

    print("[INFO] Tabix index done")
    print(
        f"[INFO] Compressed output files:\n"
        f" - {gz_tsv}\n"
        f" - {tbi_path}"
    )

    if os.path.exists(gz_tsv) and os.path.exists(tbi_path):
        os.remove(input_file)
        print(f"[INFO] Removed original uncompressed TSV: {input_file}")

    return {
        "gz_path": gz_tsv,
        "tbi_path": tbi_path,
    }