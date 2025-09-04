#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import json
from pathlib import Path
from jinja2 import Template


def main(libs: list[str], metrics: list[Path], template_path: Path, out_dir: Path):
    # Load the template
    with open(template_path) as f:
        template = Template(f.read())
    library_reads = []
    for lib, metric in zip(libs, metrics):
        # Load the metrics file
        with open(metric) as f:
            demux_metrics = json.load(f)
        passing_reads = passing_perc = error_reads = error_perc = 0.0
        for status in demux_metrics['reads']:
            if status in ['Pass', 'TooShortError']:
                passing_reads += demux_metrics['reads'][status][0]
                passing_perc += float(demux_metrics['reads'][status][1][:-1])  # Remove the '%' sign
            else:
                error_reads += demux_metrics['reads'][status][0]
                error_perc += float(demux_metrics['reads'][status][1][:-1])  # Remove the '%' sign
        library_reads.append({
            'lib_name': lib,
            'passing_reads': int(passing_reads),
            'passing_perc': f"{passing_perc:.1f}%",
            'error_reads': int(error_reads),
            'error_perc': f"{error_perc:.1f}%",
        })
    report = template.render(library_reads=library_reads)
    # Save to file
    with open("library.report.html", "w") as f:
        f.write(report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create library report')
    parser.add_argument('--libs', nargs='+', type=str, required=True, help='List of library names')
    parser.add_argument('--metrics', nargs='+', type=Path, required=True, help='List of paths to metrics json files')
    parser.add_argument('--template', type=Path, help='Library report html template')
    parser.add_argument('--outDir', type=Path, help='Directory for outputs')
    args = parser.parse_args()
    args.outDir.mkdir(parents=True, exist_ok=True)
    main(
        args.libs,
        args.metrics,
        args.template,
        args.outDir,
    )