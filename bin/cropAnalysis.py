#!/usr/bin/env python


import argparse
import collections
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


#Find any guideTab.csv files in the provided directory recursively, and save the paths in a list.
#Read in full scale RNA-seq data summary stats
def main(sample: str, guideTabCsv: Path, allcellsCsv: Path, guideMetricsCsv: Path, thresh: int, outDir: Path):
    
    allGuideTab_dict = {}
    allGuideTab_df = pd.read_csv(guideTabCsv, index_col='Cell_Barcode')
    allGuideTab_dict[sample] = allGuideTab_df

    rna = pd.read_csv(allcellsCsv, index_col=[0])
    if 'i7_alias' in rna:
       rna = rna.rename(columns={"i7_alias": "i5_alias"}) #easy merging for downstream filtering in case a legacy i7 plate was used    
    if 'I7_alias' in rna:
       rna = rna.rename(columns={"I7_alias": "i5_alias"})

    #This reads in the countGuides.py metrics list rather than calculating it in sheet, should be identical to the in-sheet generated list

    cellStats_dict = {}
    #can replace with guide list now that there isn't a mismatch in dimension between guideTab.csv and cellMetrics.csv
    cellStats_dict[sample] = pd.read_csv(guideMetricsCsv, index_col='Cell_Barcode')


    for k, df in cellStats_dict.items():
        #df = allGuideTab_dict[k]
        #print(tmp.shape)
        ax = sns.kdeplot(df['CROP_nUMI'], clip=[0,50])
        ax.set_xlabel("Total Guide's Detected (nUMI)")
        ax.set_title(sample)
        fig = ax.get_figure()
        fig.savefig(outDir / "nUMI_Distribution.png")
        plt.close()

    #Join the RNA cell stats and the CROP cell stats

    rnaGuide_dict = {}
    for k in cellStats_dict.keys():
        rnaGuide_dict[k] = rna.join(cellStats_dict[k], rsuffix='_rna', how='left')



    for k in rnaGuide_dict.keys():
        ax = plt.gca()
        rnaGuide_df = rnaGuide_dict[k]
        ax = sns.scatterplot(x=rnaGuide_df.umis+1, y=rnaGuide_df['CROP_max']+1, alpha=0.2)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel("RNA count")
        ax.set_ylabel("Guide count")
        ax.set_title(sample)
        fig = ax.get_figure()
        fig.savefig(outDir / "Guide_vs_RNA_Counts_Scatter.png")
        plt.close()

    for k, df in cellStats_dict.items():
        ax = plt.gca()
        ax = sns.scatterplot(x=df['CROP_nReads'], y=df['CROP_Saturation'], alpha=0.2, hue=(df['CROP_assignedGuide'] != 'NONE'))
        ax.set_xlabel("CROP nReads")
        ax.set_ylabel("Saturation")
        ax.set_title(sample)
        #sns.kdeplot(x=df['CROP_nReads'], y=df['CROP_Saturation'], levels=5,fill=True,alpha=0.6,cut=2)
        fig = ax.get_figure()
        fig.savefig(outDir / "Saturation_vs_nReads_Scatter.png")
        plt.close()

    for k, df in cellStats_dict.items():
        ax = plt.gca()
        ax = sns.scatterplot(x=df['CROP_nReads'], y=df['CROP_nUMI'], alpha=0.07, hue=(df['CROP_assignedGuide'] != 'NONE'))
        ax.set_xlabel("CROP nReads")
        ax.set_ylabel("CROP nUMI")
        ax.set_title(sample)
        #sns.kdeplot(x=df['CROP_nReads'], y=df['CROP_Saturation'], levels=5,fill=True,alpha=0.6,cut=2)
        fig = ax.get_figure()
        fig.savefig(outDir / "nUMI_vs_nReads_Scatter.png")
        plt.close()

    for key, df in cellStats_dict.items():
        #df = allGuideTab_dict[k]
        #print(tmp.shape)
        ax = sns.kdeplot(df['CROP_max'], clip=[0,50])
        ax.set_title(key)
        fig = ax.get_figure()
        fig.savefig(outDir / "Max_Guide_UMI_Density.png")
        plt.close()

    # Loop through the dictionary and generate subplots
    for i, (key, df) in enumerate(cellStats_dict.items()):
        ax = sns.ecdfplot(df['CROP_max'], label='Top Guides')
        sns.ecdfplot(df['CROP_nUMI'], ax=ax, label="nUMI Guides")
        ax.set_xlim(0,100)
        ax.legend()
        ax.set_title(key)
        fig = ax.get_figure()
        fig.savefig(outDir / "Proportion_of_Max_Guide.png")
        plt.close()
    # Show the plot
   


    # Loop through the dictionary and generate subplots
    for i, (key, df) in enumerate(cellStats_dict.items()):
        ax = sns.countplot(x=df.CROP_guides, label = key)
        ax.set_xlabel(f"Number of guides >= {thresh} UMIs")
        ax.set_ylabel("Cells")
        ax.set_xlim(-0.75,10)
        ax.set_title(key)
        fig = ax.get_figure()
        fig.savefig(outDir / "Cell_Count_Containing_Passing_Guides.png")
        plt.close()



    for key, df in allGuideTab_dict.items():
        ax = sns.histplot((df>=thresh).sum(0))
        ax.set_title(sample)
        ax.set_xlabel("Cells guide is found in")
        ax.set_ylabel("Guides")
        fig = ax.get_figure()
        fig.savefig(outDir / "Passing_Guide_Cell_Histogram.png")
        plt.close()


    #fil = tab[tab['sum']>=10]
    #ax = sns.kdeplot(np.log(fil['purity']+0.01))

    for key, df in cellStats_dict.items():
        ax = sns.kdeplot(df['CROP_purity'])
        ax.set_title(key)
        fig = ax.get_figure()
        fig.savefig(outDir / "Guide_Purity_Density.png")
        plt.close()


    #calculate and plot different thresholded passing densities, with a abline at the set threshold

    # Define the range of thresholds
    thresholds = range(51)

    # Initialize an empty DataFrame to store the results
    plotDF = pd.DataFrame(columns=['Threshold', 'PassingGuidesPercent', 'PassingCellsPercent', 'PassingCellsTwoGuidesPercent'])

    # Calculate percentages for each threshold
    column_percentage = []

    row_percentage = []

    row_2plus_percentage = []

    for threshold in thresholds:
        # Calculate column percentage
        column_percentage.append((allGuideTab_dict[sample] >= threshold).any(axis = 0).mean() * 100)

        # Calculate row percentage
        row_percentage.append((allGuideTab_dict[sample] >= threshold).any(axis=1).mean() * 100)

        row_2plus_percentage.append(((allGuideTab_dict[sample] >= threshold).sum(axis=1) >= 2).mean() * 100)

    # Add the results to the result DataFrame
    plotDF["Threshold"] = thresholds
    plotDF["PassingGuidesPercent"] = column_percentage
    plotDF["PassingCellsPercent"] = row_percentage
    plotDF["PassingCellsTwoGuidesPercent"] = row_2plus_percentage


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

    cellstats_df = cellStats_dict[sample]
    assigned_guides_for_count = cellstats_df["CROP_assignedGuide"].str.split(";").explode() #generate a series of all of the guide labels that have been assigned, including multiplets split into their parts
    guide_cell_counts = assigned_guides_for_count.value_counts()  #dataframe of counts of passing guides
    guide_cell_counts = guide_cell_counts.drop(labels=['NONE'], errors='ignore')
    guide_cell_counts.to_csv(outDir / "Guide_Cell_Numbers.txt", sep="\t" )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Count guide sequences matches from bcParser output')
    parser.add_argument('--sample', metavar='SAMPLE', nargs='?', const="SAMPLE", type=str, default="SAMPLE",
                        help='Sample name')
    parser.add_argument('--guideTabCsv', metavar='GUIDETAB.csv', type=Path,
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
    main(args.sample, args.guideTabCsv, args.allcellsCsv, args.guideMetricsCsv, args.thresh, args.outDir)