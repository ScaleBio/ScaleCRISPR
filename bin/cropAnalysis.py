#!/usr/bin/env python


import argparse
import collections
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.sparse import csr_matrix
from scipy.io import mmread

#Find any guideTab.csv files in the provided directory recursively, and save the paths in a list.
#Read in full scale RNA-seq data summary stats

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



def main(sample: str, umi_matrix: Path, allcellsCsv: Path, guideMetricsCsv: Path, thresh: int, outDir: Path):
    
    guide_matrix = mmread(umi_matrix)

    rna = pd.read_csv(allcellsCsv, index_col=[0])
    if 'i7_alias' in rna:
       rna = rna.rename(columns={"i7_alias": "i5_alias"}) #easy merging for downstream filtering in case a legacy i7 plate was used    
    if 'I7_alias' in rna:
       rna = rna.rename(columns={"I7_alias": "i5_alias"})

    #This reads in the countGuides.py metrics list rather than calculating it in sheet, should be identical to the in-sheet generated list

    #can replace with guide list now that there isn't a mismatch in dimension between guideTab.csv and cellMetrics.csv
    cellStats_df = pd.read_csv(guideMetricsCsv, index_col='Cell_Barcode')


    ax = sns.kdeplot(cellStats_df['CROP_nUMI'], clip=[0,50])
    ax.set_xlabel("Total Guide's Detected (nUMI)")
    ax.set_title(sample)
    fig = ax.get_figure()
    fig.savefig(outDir / "nUMI_Distribution.png")
    plt.close()

    #Join the RNA cell stats and the CROP cell stats

    rnaGuide_df = rna.join(cellStats_df, rsuffix='_rna', how='left')



    ax = plt.gca()
    ax = sns.scatterplot(x=rnaGuide_df.umis+1, y=rnaGuide_df['CROP_max']+1, alpha=0.2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("RNA count")
    ax.set_ylabel("Guide count")
    ax.set_title(sample)
    fig = ax.get_figure()
    fig.savefig(outDir / "Guide_vs_RNA_Counts_Scatter.png")
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
    fig.savefig(outDir / "Saturation_vs_nReads_Scatter.png")
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
    fig.savefig(outDir / "nUMI_vs_nReads_Scatter.png")
    plt.close()


    ax = sns.kdeplot(cellStats_df['CROP_max'], clip=[0,50])
    ax.set_title(sample)
    fig = ax.get_figure()
    fig.savefig(outDir / "Max_Guide_UMI_Density.png")
    plt.close()


    ax = sns.ecdfplot(cellStats_df['CROP_max'], label='Top Guides')
    sns.ecdfplot(cellStats_df['CROP_nUMI'], ax=ax, label="nUMI Guides")
    ax.set_xlim(0,100)
    ax.legend()
    ax.set_title(sample)
    fig = ax.get_figure()
    fig.savefig(outDir / "Proportion_of_Max_Guide.png")
    plt.close()
    # Show the plot
   



    ax = sns.countplot(x=cellStats_df['CROP_guides'], label = sample)
    ax.set_xlabel(f"Number of guides >= {thresh} UMIs")
    ax.set_ylabel("Cells")
    ax.set_xlim(-0.75,10)
    ax.set_title(sample)
    fig = ax.get_figure()
    fig.savefig(outDir / "Cell_Count_Containing_Passing_Guides.png")
    plt.close()


    #All guides will be passing in all cells at thresh = 0, and so we will not generate this plot if that's the case
    if thresh != 0:
        ax = sns.histplot((guide_matrix>=thresh).sum(1).A1)
        ax.set_title(sample)
        ax.set_xlabel("Cells guide is found in")
        ax.set_ylabel("Guides")
        fig = ax.get_figure()
        fig.savefig(outDir / "Passing_Guide_Cell_Histogram.png")
        plt.close()



    ax = sns.kdeplot(cellStats_df['CROP_purity'])
    ax.set_title(sample)
    fig = ax.get_figure()
    fig.savefig(outDir / "Guide_Purity_Density.png")
    plt.close()


    #calculate and plot different thresholded passing densities, with a abline at the set threshold

    # Define the range of thresholds
    thresholds = range(51)

    # Initialize an empty DataFrame to store the results
    plotDF = build_thresholding_assignment_df(guide_matrix, thresholds)


    ax = sns.lineplot(data=plotDF, x='Threshold', y='PassingGuidesPercent', marker='o', label='Passing Guides')
    ax = sns.lineplot(data=plotDF, x='Threshold', y='PassingCellsPercent', marker='o', label='Cells with Passing Guide', ax = ax)
    sns.lineplot(data=plotDF, x='Threshold', y='PassingCellsTwoGuidesPercent', marker='o', label='Cells with Two Passing Guides', ax = ax)
    ax.set_xlabel('Threshold')
    ax.set_ylabel('Percentage')
    ax.set_title('Passing Percentage vs. Detection Threshold')
    ax.axvline(x=thresh, color='red', linestyle='--', label=f"Threshold = {thresh}")
    ax.legend()
    fig = ax.get_figure()
    fig.savefig(outDir / "Passing_Guide_and_Cell_Percent.png")
    plotDF.to_csv(outDir / "Passing_Cell_Thresholds.txt", sep="\t", index=False)
    
    #generate count table of for each guide, how many cells had a passing value for that guide. Not differentiating doublets in this case    

    assigned_guides_for_count = cellStats_df["CROP_assignedGuide"].str.split(";").explode() #generate a series of all of the guide labels that have been assigned, including multiplets split into their parts
    guide_cell_counts = assigned_guides_for_count.value_counts()  #dataframe of counts of passing guides
    guide_cell_counts = guide_cell_counts.drop(labels=['Fail'], errors='ignore')
    guide_cell_counts.to_csv(outDir / "Guide_Cell_Numbers.txt", sep="\t" )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Count guide sequences matches from bcParser output')
    parser.add_argument('--sample', metavar='SAMPLE', nargs='?', const="SAMPLE", type=str, default="SAMPLE",
                        help='Sample name')
    parser.add_argument('--umi_matrix', metavar='MATRIX.mtx', type=Path,
                        help='countGuides.py output of passing cells x guide UMI matrix')
    parser.add_argument('--allcellsCsv', metavar='ALLCELLS.csv', type=Path,
                        help='nf-rna per-cell summary statistics output to filter based on pass status')
    parser.add_argument('--guideMetricsCsv', metavar='CELLMETRICS.csv', type=Path,
                        help='countGuides.py output of cell metadata for crop cell guide UMI matrix')
    parser.add_argument('--thresh', metavar='THRESH', nargs='?', const=3, type=int, default=3,
                        help='Theshold for the number of UMIs per guide in order to call a guide detected')
    parser.add_argument('--outDir', type=Path,
                        help='Directory for outputs')
    args = parser.parse_args()
    args.outDir.mkdir(parents=True, exist_ok=True)
    main(args.sample, args.umi_matrix, args.allcellsCsv, args.guideMetricsCsv, args.thresh, args.outDir)