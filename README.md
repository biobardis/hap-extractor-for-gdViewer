# Development Status
>
> `hap-extractor-for-gdViewer` is currently under active development.
>
> The latest development code is available in the [`dev`](../../tree/dev) branch. The `main` branch is not yet intended for general use.

# hap_extractor_for_gdViewer

`hap_extractor_for_gdViewer` is a tool for extracting haplotype information from VCF and GFF3 files 
and preparing haplotype data for visualization in [gdViewer](https://github.com/biobardis/jbrowse-plugin-gdViewer).

The program integrates genomic variants and gene annotations to identify and summarize haplotypes for selected genes or genome-wide analysis.

Gene selection must be specified with `--gene-list`. Use `--gene-list all` to analyze all genes in the GFF3 file, or provide a text file containing gene IDs to analyze only those genes.

### Requirements

* Python
* pysam

### Tested environment

The current version has been tested with the following environment:

* Python 3.14.0
* pysam 0.23.3

Install the required Python dependency using:

```bash
pip install -r requirements.txt
```

or, in a conda/mamba environment:

```bash
mamba install -c bioconda pysam
```

## Command-line options

The following command-line arguments are available:

| Argument         | Required | Default | Description |
|------------------|---|---:|---|
| `--gff`          | Yes | — | Input GFF3 annotation file. |
| `--vcf`          | Yes | — | Input phased VCF file used for haplotype analysis. The VCF should be bgzip-compressed and indexed when using `.vcf.gz`. A corresponding `.tbi` or `.csi` index file should be present in the same directory. |
| `--output-dir`   | Yes | — | Output directory. Automatically generated result files are written here. The directory is created automatically if it does not already exist. |
| `--gene-list`    | Yes | — | Gene selection. Use `all` to analyze all genes in the GFF3 file, or provide a text file containing gene IDs, one gene ID per line. |
| `--all-del`      | No | `None` | Optional pre-generated deletion TSV containing the columns `chrom`, `start`, `end`, `ref`, and `alt`. If omitted, the program automatically generates an all-DEL file from the input VCF before haplotype analysis. |
| `--report-every` | No | `100000` | Number of VCF records processed between progress reports during automatic all-DEL generation. Must be at least `1`. |
| `--mode`         | No | `CDS` | Analysis mode. Available values are `GENE` and `CDS`. |
| `--flank`        | No | `2000` | Length of the flanking region, in base pairs, used in `GENE` mode. Must be `0` or greater. |
| `--min-sample-count` | No | `2` | Minimum number of occurrences required for a haplotype to be retained as a `HAP` row. Must be at least `1`. |
| `--all-del-output`     | No | Automatically generated | Optional output path for the automatically generated all-DEL TSV when `--all-del` is not supplied. |
| `--hap-output`   | No | Automatically generated | Optional filename or path for the main output TSV. If omitted, the filename is generated automatically. |
| `--summary-output` | No | Automatically generated | Optional filename or path for the summary TSV. If omitted, the filename is generated automatically. |
| `--debug`        | No | Disabled | Enable debug mode. Only the first `--debug-n-genes` genes are processed. |
| `--debug-n-genes` | No | `100` | Number of genes processed when `--debug` is enabled. Must be at least `1`. |

## Usage

Analyze selected genes using a gene-list file:

```bash
python hap_extract.py \
  --gff annotation.gff3 \
  --vcf variants.vcf.gz \
  --gene-list genes.txt \
  --output-dir results
```

The gene-list file should contain one gene ID per line:

```text
Gene001
Gene002
Gene003
```

To analyze all genes in the GFF3 file, use:

```bash
python hap_extract.py \
  --gff annotation.gff3 \
  --vcf variants.vcf.gz \
  --gene-list all \
  --output-dir results
```

Show all available options with:

```bash
python hap_extract.py -h
```

## Output

If an output file with the same name already exists, a numeric suffix is added automatically.

For output files larger than 100 MB, the program automatically generates bgzip-compressed and Tabix-indexed files.
