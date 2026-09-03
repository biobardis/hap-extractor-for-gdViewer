"""Program entry point."""

import time

from .cli import parse_args
from .config import AnalysisConfig


def main(argv=None):
    args = parse_args(argv)
    config = AnalysisConfig.from_namespace(args)

    # Import here so `-h` remains usable even before pysam is installed.
    from .pipeline import run_pipeline

    start = time.time()
    result = run_pipeline(config)
    end = time.time()

    print(f"Total runtime: {end - start:.2f} seconds")
    return result