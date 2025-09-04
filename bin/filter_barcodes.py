#!/usr/bin/env python


import argparse
from pathlib import Path
import polars as pl

from scipy.sparse import csr_array, hstack
from scipy.io import mmread, mmwrite


def filter_sparse_matrix(raw_mtx, filtered_indices, cell_stats):
    """
    Rebuild sparse matrix to contain only the cell barcodes that were called a cell in the allCells.csv from RNA.
    Cell barcodes can have passed in RNA that are not present in the raw matrix generated from the enriched library,
    and so we cannot do a simple subsetting
    Args:
        raw_mtx: raw umi or read csr sparse matrix in guides x cells orientation
        filtered_indices: passing barcodes from RNA allCells.csv
        cell_stats: cell metrics file to provide raw barcodes list
    Returns:
        filtered_mtx:
            counts matrix in guides x cells orientation where cells are only barcodes that were called a cell
            in ScaleRNA, with empty columns for those that were not detected in enriched library
    """
    # get indices in raw matrix of barcodes that were called a cell in RNA
    indices = (
        filtered_indices
        .with_row_index("filtered_index")
        .join(
            cell_stats.select(["Cell_Barcode"]).with_row_index(), left_on="cell_id", right_on="Cell_Barcode", how="left"
        )
    )
    idx = indices.filter(pl.col("index").is_not_null())["index"].to_numpy()
    filtered_mtx = raw_mtx[:, idx]

    # Add sparse columns for barcodes that were not detected in the enriched library
    missing_idx = indices.filter(pl.col("index").is_null())["filtered_index"].to_numpy()
    for idx in missing_idx:
        filtered_mtx = hstack([filtered_mtx[:, :idx], csr_array((filtered_mtx.shape[0], 1)), filtered_mtx[:, idx:]])

    return filtered_mtx

def get_crispr_to_rna_mapping(fname: Path) -> dict[str, str]:
    """Read in file for mapping CRISPR to RNA PCR barcodes"""
    with open(fname, "r") as f:
        return {line.split("\t")[0]: line.split("\t")[1].strip() for line in f}


def main(
    umi_matrix: Path,
    reads_matrix: Path,
    barcodes_tsv: Path,
    cell_metrics: Path,
    out_dir: Path,
    thresh: int,
    crispr_to_rna_mapping: Path,
):
    mapping = get_crispr_to_rna_mapping(crispr_to_rna_mapping) if crispr_to_rna_mapping else {}
    umi_matrix = mmread(umi_matrix).tocsr()
    reads_matrix = mmread(reads_matrix).tocsr()
    # filter cellMetrics to only those in the RNA data that match the passing criteria in the RNA data
    rna = pl.read_csv(barcodes_tsv, has_header=False, new_columns=["cell_id"])
    
    cell_stats_df = pl.read_csv(cell_metrics)
    # replace CRISPR PCR barcodes with RNA PCR barcodes if mapping is provided
    cell_stats_df = cell_stats_df.with_columns(
        pl.col("Cell_Barcode").str.split(by="+").list.eval(pl.element().replace(mapping)).list.join("+").alias("Cell_Barcode")
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    # write raw barcodes
    with open(out_dir / "barcodes.tsv", "w") as f:
        cell_stats_df.select(["Cell_Barcode"]).write_csv(f, separator="\t", include_header=False)

    guide_matrix_filtered = filter_sparse_matrix(umi_matrix, rna, cell_stats_df)
    with open(out_dir / "filtered.umi.matrix.mtx", "wb") as f:
        mmwrite(f, guide_matrix_filtered)
    with open(out_dir / "filtered.barcodes.tsv", "w") as f:
        rna.write_csv(f, separator="\t", include_header=False)

    reads_matrix_filtered = filter_sparse_matrix(reads_matrix, rna, cell_stats_df)
    with open(out_dir / "filtered.reads.matrix.mtx", "wb") as f:
        mmwrite(f, reads_matrix_filtered)

    # Create merged allCells with all rows from RNA and matching rows from ENRICH if present
    cell_stats_passing = rna.join(
        cell_stats_df, how="left", left_on="cell_id", right_on="Cell_Barcode", suffix="_ENRICH"
    ).fill_null(strategy="zero")
    cell_stats_passing.with_columns(
        #Set all non-assigned guide cells to 'Fail'
        pl.when(pl.col("CROP_max") < thresh).then(pl.lit("Fail")).otherwise(pl.col("CROP_assignedGuide")).alias("CROP_assignedGuide")
    )
    cell_stats_passing.write_csv(out_dir / "cellMetrics.filtered.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate filtered guide umi and read counts matrix and cellMetrics \
            that share the same cell barcodes as the passing RNA analysis per sample"
    )
    parser.add_argument(
        "--umiMatrix", type=Path, help="countGuides output of guide x cells raw UMI matrix"
    )
    parser.add_argument(
        "--readsMatrix", type=Path, help="countGuides output of guide x cells raw reads matrix"
    )
    parser.add_argument(
        "--barcodes",
        type=Path,
        help="Path to barcodes.tsv.gz file containing passing cell barcodes",
    )
    parser.add_argument(
        "--cellMetrics",
        type=Path,
        help="countGuides output of cell metadata for cell guide UMI matrix",
    )
    parser.add_argument(
        "--crisprToRnaMapping",
        type=Path,
        help="Path to crispr to rna mapping file",
    )
    parser.add_argument('--thresh', metavar='THRESH', nargs='?', const=3, type=int, default=3,
                        help='Theshold for the number of UMIs per guide in order to call a guide detected')
    parser.add_argument("--outDir", type=Path, default=".", help="Directory for outputs")
    args = parser.parse_args()
    args.outDir.mkdir(parents=True, exist_ok=True)
    main(
        args.umiMatrix,
        args.readsMatrix,
        args.barcodes,
        args.cellMetrics,
        args.outDir,
        args.thresh,
        args.crisprToRnaMapping,
    )

