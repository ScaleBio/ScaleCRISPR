#!/usr/bin/env python
"""
Merge raw umi count matrices
"""
import argparse
from pathlib import Path
from scipy.io import mmread, mmwrite
from scipy.sparse import hstack


def merge_mtx(output_file: str, mtx_files: list[Path]):
    """
    Merge raw umi matrices

    Sparse coordinate matrices are horizontally concatenated, rows are features and columns are cells.

    Args:
        output_file: Output matrix file name
        mtx_files: List of sparse coordinate mtx files to merge
    """
    raw_mtx = [mmread(file) for file in mtx_files]
    with open(output_file, "wb") as f:
        mmwrite(f, hstack(raw_mtx))


def main():
    parser = argparse.ArgumentParser("Merge individual umi count matrices for sample")
    parser.add_argument("mtx", nargs="+", type=Path, help="sparse coordinate mtx files to merge")
    parser.add_argument("--outputFile", type=Path, required=True, help="Output file name")
    args = parser.parse_args()

    merge_mtx(output_file=args.outputFile, mtx_files=args.mtx)


if __name__ == "__main__":
    main()
