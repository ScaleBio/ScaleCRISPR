#!/usr/bin/env python
"""
Concatenate columns (cells) for gene expression matrices in Matrix Market format, into single output file
"""

import duckdb
import argparse
from pathlib import Path


def concat_mtx(mtx_files: list[Path], output_file: Path, con: duckdb.DuckDBPyConnection) -> None:
    """
    Concatenate sparse matrices in Matrix Market format

    Args:
        mtx_files: List of filtered mtx files to concatenate
        output_file: Output path
        con: duckdb connection
    """
    nnz = 0
    row_size = []
    col_size = []
    for idx, mtx_file in enumerate(mtx_files):
        # Open mtx file to get dimensions and check if it represents an empty mtx
        with open(mtx_file, "rt") as f:
            # https://math.nist.gov/MatrixMarket/formats.html
            f.readline()  # %%MatrixMarket matrix coordinate real general
            f.readline()  # %
            dims = f.readline().split()  # 36910 185326 0000262179120
            nnz += int(dims[2])
            row_size.append(int(dims[0]))
            col_size.append(int(dims[1]))
        if col_size[idx] == 0:
            # Skip matrices with no cells
            continue
        con.sql(
            f"""
        INSERT INTO mtx (
            SELECT gene, barcode + {sum(col_size[:idx])} as barcode, count
            FROM read_csv(
                '{mtx_file}',
                delim=' ',
                skip={3},
                columns = {{
                    'gene': 'UINTEGER',
                    'barcode': 'UINTEGER',
                    'count': 'FLOAT'
                }})
        );
        """
        )
    if len(set(row_size)) != 1:
        raise AssertionError("Row dimensions must match to concatenate Matrix Market files by column")
    con.sql(
        f"""
    COPY mtx
    TO '{output_file}' (
        FORMAT csv,
        DELIMITER ' ',
        HEADER false,
        PREFIX '%%MatrixMarket matrix coordinate real general\n%\n{row_size[0]} {sum(col_size)} {str(nnz).zfill(13)}\n',
        SUFFIX '\n'
    );
    """
    )


def main():
    parser = argparse.ArgumentParser("Combine Matrix Market files into single output file, stack columns")
    parser.add_argument("mtxFiles", nargs="+", type=Path, help="Filtered matrix files to merge")
    parser.add_argument("--id", type=str, help="Sample ID")
    parser.add_argument("--outputFile", required=True, type=Path, help="Output path for merged file")
    parser.add_argument("--threads", type=int, required=False, default=2, help="Number of threads for duckdb")
    parser.add_argument("--memory", type=str, required=False, default="7 GB", help="Memory allocated to task")
    args = parser.parse_args()

    mem_limit, mem_unit = args.memory.split()
    # allocate duckdb memory based on number of threads
    mem_limit = f"{float(mem_limit) / (args.threads + 1):.1f}{mem_unit}"
    # create persistent duckdb connection
    con = duckdb.connect("tables.db")
    con.sql(
        f"""
    SET threads TO {args.threads};
    SET memory_limit TO '{mem_limit}';
    CREATE OR REPLACE TABLE mtx (
        gene UINTEGER,
        barcode UINTEGER,
        count FLOAT
    );
    """
    )
    concat_mtx(mtx_files=args.mtxFiles, output_file=args.outputFile, con=con)
    con.close()


if __name__ == "__main__":
    main()
