#!/usr/bin/env python
"""
Create sample-level statistics from cell metrics
"""

import argparse
from pathlib import Path
from dataclasses import dataclass
import duckdb


@dataclass
class Metrics:
    """Overall read statistics for the dataset"""
    reads: int = 0 # Input reads
    barcodeError: int = 0 # Reads without a valid cell barcode
    noGuide: int = 0 # Reads without a CRISPR Guide sequence match
    correct: int = 0 # Usable reads
    meanReadsPerCell: int = 0 # number of reads per cell. Should be multiplied by correct for ~usable amount
    nUMIPerCell: int = 0 # number of CROP UMI molecules per cell
    meanUMIPerGuide: int = 0 # mean number of UMI associated with each guide across the cells
    passingPercent: int = 0 # percent of cells that had at least thresh Guide UMI's detected

    def print(self, outFn: Path):
        with open(outFn, 'w') as out:
            print("Reads", f"{self.reads:}", sep=',', file=out)
            print("NoGuide", f"{self.noGuide/self.reads:.1%}", sep=',', file=out)
            print("Correct", f"{self.correct/self.reads:.1%}", sep=',', file=out)
            print("ReadsPerCell", f"{self.meanReadsPerCell:}", sep=',', file=out)
            print("nUMIPerCell", f"{self.nUMIPerCell:}", sep=',', file=out)
            print("meanUMIPerGuide",f"{self.meanUMIPerGuide:}", sep=',', file=out)
            print("passingPercent", f"{self.passingPercent:.1%}", sep=',', file=out)


def calculate_sample_stats(
    cell_metrics: list[Path],
    filtered_cell_metrics: list[Path],
    num_guides: int,
    thresh: int,
) -> Metrics:
    """
    Calculate sample-level statistics from the cellMetrics csv files

    Args:
        cellMetrics: List of paths to the raw cellMetrics csv files
        cellMetricsFiltered: List of paths to the cellMetrics csv files aligned to passing RNA cells
        num_guides: Number of guides for this sample
        thresh: Threshold for the number of UMIs per guide in order to call a guide detected

    Returns:
        Sample-level statistics
    """
    metrics = Metrics()

    cell_stats = (
        duckdb.sql(
            f"""
        SELECT
            sum(CROP_nReads) AS reads,
            sum(CROP_noGuide * CROP_nReads) AS no_guide,
        FROM read_csv({[str(x) for x in cell_metrics]})
        """
        )
        .df()
        .to_dict(orient="index")[0]
    )

    cell_stats_passing = (
        duckdb.sql(
            f"""
        SELECT
            mean(CROP_nReads) AS mean_reads_per_cell,
            median(CROP_nUMI) AS n_umi_per_cell,
            mean(CAST(CROP_max >= {thresh} AS INTEGER)) AS passing_percent,
            sum(CROP_nUMI) AS total_umi,
        FROM read_csv({[str(x) for x in filtered_cell_metrics]})
        """
        )
        .df()
        .to_dict(orient="index")[0]
    )

    metrics.reads = round(cell_stats['reads'])
    metrics.noGuide = cell_stats['no_guide']
    metrics.correct = metrics.reads - metrics.noGuide
    metrics.meanReadsPerCell = round(cell_stats_passing['mean_reads_per_cell'], 1)
    metrics.nUMIPerCell = round(cell_stats_passing['n_umi_per_cell'], 1)
    metrics.passingPercent = cell_stats_passing['passing_percent']
    metrics.meanUMIPerGuide = round(cell_stats_passing['total_umi'] / num_guides, 1)

    return metrics


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--cellMetrics", nargs="+", type=Path, required=True, help="Path to the raw cellMetrics csv files"
    )
    parser.add_argument(
        "--cellMetricsFiltered", nargs="+", type=Path, required=True, help="Path to the cellMetrics csv files aligned to passing RNA cells"
    )
    parser.add_argument("--sample", type=str, required=True, help="Unique string to identify this sample.")
    parser.add_argument("--numGuides", type=int, required=False, default=1, help="Number of guides for this sample")
    parser.add_argument('--thresh', metavar='THRESH', nargs='?', const=3, type=int, default=3,
                        help='Theshold for the number of UMIs per guide in order to call a guide detected')
    parser.add_argument("--threads", type=int, required=False, default=1, help="Number of threads for duckdb")
    parser.add_argument("--memory", type=str, required=False, default="8 GB", help="Memory allocated to task")

    args = parser.parse_args()

    mem_limit, mem_unit = args.memory.split()
    # allocate duckdb memory based on number of threads
    mem_limit = f"{float(mem_limit) / (args.threads + 1):.1f}{mem_unit}"
    duckdb.sql(
        f"""
    SET threads TO {args.threads};
    SET memory_limit TO '{mem_limit}';
    """
    )

    sample_stats = calculate_sample_stats(
        cell_metrics=args.cellMetrics,
        filtered_cell_metrics=args.cellMetricsFiltered,
        num_guides=args.numGuides,
        thresh=args.thresh,
    )

    metrics_dir = Path(".")
    sample_stats.print(metrics_dir / "metrics.csv")


if __name__ == "__main__":
    main()
