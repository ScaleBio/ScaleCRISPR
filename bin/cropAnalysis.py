#!/usr/bin/env python


import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.sparse import csr_matrix
from scipy.io import mmread
from jinja2 import Template
import io
import base64
import datetime

#Find any guideTab.csv files in the provided directory recursively, and save the paths in a list.
#Read in full scale RNA-seq data summary stats

def plot_to_base64():
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode('utf-8')
    return f"data:image/png;base64,{encoded}"


def build_thresholding_assignment_df(sparse_matrix, thresholds):

    """Create the dataframe that captures the percentage of passing cells and guides per threshold for plotting"""
    # Initialize an empty DataFrame to store the results
    plotDF = pd.DataFrame(columns=['Threshold', 'PassingGuidesPercent', 'PassingCellsPercent', 'PassingCellsTwoGuidesPercent'])

    # Calculate percentages for each threshold
    column_percentage = []

    row_percentage = []

    column_2plus_percentage = []

    num_columns = sparse_matrix.shape[1]
    num_rows = sparse_matrix.shape[0]

    for threshold in thresholds:
        # Calculate column percentage
        if threshold == 0:
            column_percentage.append(100)
            row_percentage.append(100)
            column_2plus_percentage.append(100)
        else:

            column_percentage.append(((sparse_matrix >= threshold).sum(axis = 0).astype(bool).sum() / num_columns) * 100)
            column_2plus_percentage.append((((sparse_matrix >= threshold).sum(axis=0) >= 2).sum() / num_columns) * 100)

            # Calculate row percentage
            row_percentage.append(((sparse_matrix >= threshold).sum(axis=1).astype(bool).sum() / num_rows) * 100)
  

    # Add the results to the result DataFrame
    plotDF["Threshold"] = thresholds
    plotDF["PassingGuidesPercent"] = row_percentage
    plotDF["PassingCellsPercent"] = column_percentage
    plotDF["PassingCellsTwoGuidesPercent"] = column_2plus_percentage

    return plotDF


def main(sample: str, umi_matrix: Path, guideMetricsCsv: Path, thresh: int, template_path: Path, outDir: Path, lib_json_name: str, barcodes: str, workflow_version: str):
    
    guide_matrix = mmread(umi_matrix)
    plots = []

    #This reads in the countGuides.py metrics list rather than calculating it in sheet, should be identical to the in-sheet generated list

    #can replace with guide list now that there isn't a mismatch in dimension between guideTab.csv and cellMetrics.csv
    cellStats_df = pd.read_csv(guideMetricsCsv, index_col='cell_id')

    #calculate and plot different thresholded passing densities, with a abline at the set threshold

    # Define the range of thresholds
    thresholds = range(51)

    # Initialize an empty DataFrame to store the results
    plotDF = build_thresholding_assignment_df(guide_matrix, thresholds)

    # Passing percentage vs detection threshold plot
    ax = sns.lineplot(data=plotDF, x='Threshold', y='PassingGuidesPercent', marker='o', label='Passing Guides')
    ax = sns.lineplot(data=plotDF, x='Threshold', y='PassingCellsPercent', marker='o', label='Cells with Passing Guide', ax = ax)
    sns.lineplot(data=plotDF, x='Threshold', y='PassingCellsTwoGuidesPercent', marker='o', label='Cells with Two Passing Guides', ax = ax)
    ax.set_xlabel('Threshold')
    ax.set_ylabel('Percentage')
    ax.set_title('Passing Percentage vs. Detection Threshold')
    ax.axvline(x=thresh, color='red', linestyle='--', label=f"Threshold = {thresh}")
    ax.legend()
    fig = ax.get_figure()
    fig.savefig(outDir / "plots" / f"{sample}_Passing_Guide_and_Cell_Percent.png")
    plots.append({"title": sample, "src": plot_to_base64()})
    plt.close()
    plotDF.to_csv(outDir / "csv" / f"{sample}_Passing_Cell_Thresholds.txt", sep="\t", index=False)
    
    # Density of CROP_purity plot
    ax = sns.kdeplot(cellStats_df['CROP_purity'])
    ax.set_title(sample)
    fig = ax.get_figure()
    fig.savefig(outDir / "plots" / f"{sample}_Guide_Purity_Density.png")
    plots.append({"title": sample, "src": plot_to_base64()})
    plt.close()
    
    # Histogram of cells by number of guides >= 3 UMIs plot
    ax = sns.countplot(x=cellStats_df['CROP_guides'], label = sample)
    ax.set_xlabel(f"Number of guides >= {thresh} UMIs")
    ax.set_ylabel("Cells")
    ax.set_xlim(-0.75,10)
    ax.set_title(sample)
    fig = ax.get_figure()
    fig.savefig(outDir / "plots" / f"{sample}_Cell_Count_Containing_Passing_Guides.png")
    plots.append({"title": sample, "src": plot_to_base64()})
    plt.close()


    ax = sns.kdeplot(cellStats_df['CROP_nUMI'], clip=[0,50])
    ax.set_xlabel("Total Guides Detected (nUMI)")
    ax.set_title(sample)
    fig = ax.get_figure()
    fig.savefig(outDir / "plots" / f"{sample}_nUMI_Distribution.png")
    plots.append({"title": sample, "src": plot_to_base64()})
    plt.close()

    ax = plt.gca()
    ax = sns.scatterplot(x=cellStats_df['CROP_nReads'], y=cellStats_df['CROP_Saturation'], alpha=0.2, hue=(cellStats_df['CROP_assignedGuide'] != 'Fail'))
    ax.set_xlabel("CROP nReads")
    ax.set_ylabel("Saturation")
    leg = plt.legend()
    for lh in leg.legend_handles:
            lh.set_alpha(1)
    ax.set_title(sample)
    fig = ax.get_figure()
    fig.savefig(outDir / "plots" / f"{sample}_Saturation_vs_nReads_Scatter.png")
    plots.append({"title": sample, "src": plot_to_base64()})
    plt.close()

    ax = plt.gca()
    ax = sns.scatterplot(x=cellStats_df['CROP_nReads'], y=cellStats_df['CROP_nUMI'], alpha=0.07, hue=(cellStats_df['CROP_assignedGuide'] != 'Fail'))
    ax.set_xlabel("CROP nReads")
    leg = plt.legend()
    for lh in leg.legend_handles:
            lh.set_alpha(1)
    ax.set_ylabel("CROP nUMI")
    ax.set_title(sample)
    fig = ax.get_figure()
    fig.savefig(outDir / "plots" / f"{sample}_nUMI_vs_nReads_Scatter.png")
    plots.append({"title": sample, "src": plot_to_base64()})
    plt.close()

    ax = sns.kdeplot(cellStats_df['CROP_max'], clip=[0,50])
    ax.set_title(sample)
    fig = ax.get_figure()
    fig.savefig(outDir / "plots" / f"{sample}_Max_Guide_UMI_Density.png")
    plots.append({"title": sample, "src": plot_to_base64()})
    plt.close()

    ax = sns.ecdfplot(cellStats_df['CROP_max'], label='Top Guides')
    sns.ecdfplot(cellStats_df['CROP_nUMI'], ax=ax, label="nUMI Guides")
    ax.set_xlim(0,100)
    ax.legend()
    ax.set_title(sample)
    fig = ax.get_figure()
    fig.savefig(outDir / "plots" / f"{sample}_Proportion_of_Max_Guide.png")
    plots.append({"title": sample, "src": plot_to_base64()})
    plt.close()

    #All guides will be passing in all cells at thresh = 0, and so we will not generate this plot if that's the case
    if thresh != 0:
        ax = sns.histplot((guide_matrix>=thresh).sum(1).A1)
        ax.set_title(sample)
        ax.set_xlabel("Cells guide is found in")
        ax.set_ylabel("Guides")
        fig = ax.get_figure()
        fig.savefig(outDir / "plots" / f"{sample}_Passing_Guide_Cell_Histogram.png")
        plots.append({"title": sample, "src": plot_to_base64()})
        plt.close()

    #generate count table of for each guide, how many cells had a passing value for that guide. Not differentiating doublets in this case    

    assigned_guides_for_count = cellStats_df["CROP_assignedGuide"].str.split(";").explode() #generate a series of all of the guide labels that have been assigned, including multiplets split into their parts
    guide_cell_counts = assigned_guides_for_count.value_counts()  #dataframe of counts of passing guides
    guide_cell_counts = guide_cell_counts.drop(labels=['Fail'], errors='ignore')
    guide_cell_counts.to_csv(outDir / "csv" / f"{sample}_Guide_Cell_Numbers.txt", sep="\t" )

    # Generate a HTML report from the template
    sample_info = [
         ("Sample Name", sample.split(".")[0]),
         ("Library Name", sample.split(".")[1]),
         ("Sample Barcodes", barcodes),
         ("Workflow Version", workflow_version),
         ("Analysis Date", datetime.datetime.now().strftime("%Y-%m-%d %H:%M")),
         ("Library Structure", lib_json_name)
    ]
    metrics = pd.read_csv("metrics.csv", header=None)
    read_metrics_names = ["Reads", "NoGuide", "Correct"]
    read_metrics = list(metrics[metrics.iloc[:, 0].isin(read_metrics_names)].itertuples(index=False, name=None))
    cell_metrics = list(metrics[~metrics.iloc[:, 0].isin(read_metrics_names)].itertuples(index=False, name=None))
    # Load the template
    with open(template_path) as f:
        template = Template(f.read())
    report = template.render(sample_info=sample_info, read_metrics=read_metrics, cell_metrics=cell_metrics, plots=plots)
    # Save to file
    with open(f"{sample}.report.html", "w") as f:
        f.write(report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create sample report')
    parser.add_argument('--sample', metavar='SAMPLE', nargs='?', const="SAMPLE", type=str, default="SAMPLE",
                        help='Sample name')
    parser.add_argument('--umi_matrix', metavar='MATRIX.mtx', type=Path,
                        help='countGuides.py output of passing cells x guide UMI matrix')
    parser.add_argument('--guideMetricsCsv', metavar='CELLMETRICS.csv', type=Path,
                        help='countGuides.py output of cell metadata for crop cell guide UMI matrix')
    parser.add_argument('--thresh', metavar='THRESH', nargs='?', const=3, type=int, default=3,
                        help='Theshold for the number of UMIs per guide in order to call a guide detected')
    parser.add_argument('--template', type=Path, help='QC report html template')
    parser.add_argument('--outDir', type=Path, help='Directory for outputs')
    parser.add_argument('--libJsonName', type=str, help='Name of libJson for report')
    parser.add_argument('--barcodesRange', type=str, help='RT wells used')
    parser.add_argument('--workflowVersion', type=str, help='Workflow version')
    args = parser.parse_args()
    args.outDir.mkdir(parents=True, exist_ok=True)
    (args.outDir / "csv").mkdir(parents=True, exist_ok=True)
    (args.outDir / "plots").mkdir(parents=True, exist_ok=True)
    main(
        args.sample,
        args.umi_matrix,
        args.guideMetricsCsv,
        args.thresh,
        args.template,
        args.outDir,
        args.libJsonName,
        args.barcodesRange,
        args.workflowVersion,
    )