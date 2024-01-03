#!/usr/bin/env python


import argparse
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Tuple
from collections import defaultdict
import sys
import csv
import collections
import numpy as np
import pandas as pd

# Positive matches in bcParser output
MatchStatus = ("exact", "corrected")

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
            #print("BarcodeError", f"{self.barcodeError/self.reads:.1%}", sep=',', file=out)
            print("NoGuide", f"{self.noGuide/self.reads:.1%}", sep=',', file=out)
            print("Correct", f"{self.correct/self.reads:.1%}", sep=',', file=out)
            print("ReadsPerCell", f"{self.meanReadsPerCell:}", sep=',', file=out)
            print("nUMIPerCell", f"{self.nUMIPerCell:}", sep=',', file=out)
            print("meanUMIPerGuide",f"{self.meanUMIPerGuide:}", sep=',', file=out)
            print("passingPercent", f"{self.passingPercent:.1%}", sep=',', file=out)

class cellReadCounts:
    """Track reads counts for cell (e.g cell X read combination)"""
    counts: Dict[Tuple[str, str, str], int] # Cell Sequence -> ReadCount
    def __init__(self):
        self.counts = {}
    
    def CROP_nReads(self, read_thres: int):
        return sum(1 for i in self.counts.values() if i >= read_thres)
    
    def addCount(self, cellSeq: Tuple[str, str, str]):
        self.counts[cellSeq] = self.counts.get(cellSeq, 0) + 1

class UmiReadCounts:
    """Track reads counts for each UMI within one target (e.g cell X guide combination)"""
    counts: Dict[str, int] # UMI-sequence -> ReadCount
    def __init__(self):
        self.counts = {}
    
    def nUmis(self, read_thres: int):
        return sum(1 for i in self.counts.values() if i >= read_thres)
    
    def addRead(self, umiSeq: str):
        self.counts[umiSeq] = self.counts.get(umiSeq, 0) + 1

class GuideReadCounts:
    """UMI counts for each guide in a cell"""
    guides: Dict[str, UmiReadCounts] # GuideSequence -> UMI-Counts
    def __init__(self):
        self.guides = {}
    
    def umiCounts(self, read_thres: int):
        res = {}
        for guide in self.guides:
            res[guide] = self.guides[guide].nUmis(read_thres)
        return res
    
    def guideReadCounts(self, umi: str):
        res = {}
        for guide in self.guides:
            res[guide] = self.guides[guide].counts.get(umi,0)
        return res

    def addRead(self, guide: str, umi: str):
        if not guide in self.guides:
            self.guides[guide] = UmiReadCounts()
        self.guides[guide].addRead(umi)


def countGuideReads(barcodesCsv: Path):
    """Count reads for each cell X guide X UMI combo"""
    # BC -> (GUIDE,UMI) -> ReadCount
    # BC is a tuple of three strings for the three cell-barcode levels (RT, Ligation, PCR)
    cellGuideReads: Dict[Tuple[str, str, str], GuideReadCounts] = collections.defaultdict(GuideReadCounts)
    cellReads = cellReadCounts()
    metrics = Metrics()
    for row in csv.DictReader(open(barcodesCsv),delimiter="\t"):
        metrics.reads += 1
        if not (row['lig_info'].split(',')[0] in MatchStatus and
        row['rt_info'].split(',')[0] in MatchStatus and
        row['pcr_info'].split(',')[0] in MatchStatus):
                metrics.barcodeError += 1
                continue
        cell = (row['lig'],row['rt'], row['pcr'])
        cellReads.addCount(cell) #Count a read for the cell, having had a correct match to the barcode info
        if not row['guide_info'].split(',')[0] in MatchStatus:
            metrics.noGuide += 1
            cellGuideReads[cell].addRead('noGuide',row['Read']) #Count reads that do not have a guide sequence but did have a correct cellular barcode
            continue
        metrics.correct += 1
        cellGuideReads[cell].addRead(row['guide'],row['umi'])
    return cellGuideReads, metrics, cellReads

def guideAssignment(guideTab,guideThresh):
    mxs = guideTab.ge(guideThresh, axis=0) #any column that has a value >= to the threshold, for each cell
    # join the column names of the passing values of each row into a single string separated by a ";" character
    return mxs.dot(mxs.columns + ';').str.rstrip(';')


def cellBarcodeConcat(referencePath,inputDF):
    ligation_df = pd.read_csv(f'{referencePath}/3lvlRNA_lig.txt', sep="\t", names=['Barcode','Alias'])
    ligation_dict = dict(zip(ligation_df.Alias,ligation_df.Barcode))
    rt_df = pd.read_csv(f'{referencePath}/3lvlRNA_rt.txt', sep="\t", names=['Barcode','Alias'])
    rt_dict = dict(zip(rt_df.Alias,rt_df.Barcode))
    pcr_df = pd.read_csv(f'{referencePath}/3lvlRNA_pcr.txt', sep="\t", names=['Barcode','Alias'])
    pcr_dict = dict(zip(pcr_df.Alias,pcr_df.Barcode))
    alias_df = pd.DataFrame(inputDF.index.tolist(), columns = ['Ligation_alias', 'RT_alias', 'i5_alias'])
    alias_df['Ligation_alias'] = alias_df['Ligation_alias'].replace(ligation_dict)
    alias_df['RT_alias'] = alias_df['RT_alias'].replace(rt_dict)
    alias_df['i5_alias'] = alias_df['i5_alias'].replace(pcr_dict)
    alias_df['Cell_Barcode'] = alias_df['Ligation_alias'].astype(str) + alias_df['RT_alias'].astype(str) + alias_df['i5_alias'].astype(str)
    inputDF.insert(0,"Cell_Barcode",alias_df['Cell_Barcode'].values)


def main(barcodesCsv: Path, allcellsCsv: Path, recountGuideTab: Path, recountCellMetrics: Path, recountGuideReads: Path, thresh: int, references: Path, outDir: Path):
    if (recountGuideTab is not None) and (recountCellMetrics is not None):
        cellGuideTab = pd.read_csv(recountGuideTab, index_col=[0,1,2])
        cellMetrics = pd.read_csv(recountCellMetrics, index_col=[0,1,2])
        cellGuideReads = pd.read_csv(recountGuideReads, index_col=[0,1,2])
        metrics = Metrics()
        metrics.reads = 1
    else:
        cellGuideReads, metrics, cellReads = countGuideReads(barcodesCsv)
        # Create  a table with UMI counts for each cell X guide combo
        # cellGuideCounts: Dict[]= {}
        cellMetrics = pd.DataFrame.from_dict(cellReads.counts, orient='index', columns=['CROP_nReads'])
        cellMetrics.index = pd.MultiIndex.from_tuples(cellMetrics.index, names=['Ligation_alias', 'RT_alias', 'i5_alias'])
        
        cellGuideReadCounts = defaultdict()
        cellGuideCounts = defaultdict()
        for cell in cellGuideReads:
            cellGuideCounts[cell] = {}
            cellGuideReadCounts[cell] = {}
            for guide in cellGuideReads[cell].guides:
                cellGuideCounts[cell][guide] = cellGuideReads[cell].guides[guide].nUmis(1)
                cellGuideReadCounts[cell][guide] = sum(cellGuideReads[cell].guides[guide].counts.values())
        
        cellGuideReads = pd.DataFrame.from_dict(cellGuideReadCounts).transpose().fillna(0)
        cellGuideReads.index.names = ['Ligation_alias', 'RT_alias', 'i5_alias']

        cellGuideTab = pd.DataFrame.from_dict(cellGuideCounts).transpose().fillna(0)
        cellGuideTab.index.names = ['Ligation_alias', 'RT_alias', 'i5_alias']

        #Continue updating cellMetrics df

        cellMetrics['CROP_noGuide'] = cellGuideTab['noGuide']/cellMetrics['CROP_nReads'] #Percentage of reads with no guide detected
        cellGuideTab = cellGuideTab.drop(['noGuide'], axis=1) # remove the noGuide column for downstream calcualtions and ultimate writing out of guide matrix
        cellMetrics = cellMetrics.filter(items = cellGuideTab.index, axis = 0) # reconcile any cells that might have discrepancy between the two after this filtering.
        cellMetrics['CROP_nUMI'] = (cellGuideTab.sum(1)) # number of unique guide UMI detected 
        cellMetrics['CROP_nGuide'] = (cellGuideTab>0).sum(1) # number of unique guide's detected, like unique genes detected
        cellMetrics['CROP_max'] = cellGuideTab.max(1) # max number of UMI detected across the guides within a cell
        cellMetrics['CROP_second'] = cellGuideTab.apply(lambda row: row.nlargest(2).values[-1],axis=1) # second highest number of UMI of any guide within a cell
        cellMetrics['CROP_purity'] = (cellMetrics['CROP_max'] / cellMetrics['CROP_nUMI']) 
        cellMetrics['CROP_topTwo'] = ((cellMetrics['CROP_max']+cellMetrics['CROP_second']) / cellMetrics['CROP_nUMI'])
        cellMetrics['CROP_minorFrac'] = cellMetrics['CROP_second'] / cellMetrics['CROP_max']
        cellMetrics['CROP_guides'] = (cellGuideTab>=thresh).sum(1) # how many guides per cell passed a minimum detection threshold
        cellMetrics['CROP_Saturation'] = 1 - (cellMetrics.CROP_nUMI / (cellMetrics.CROP_nReads - (cellMetrics.CROP_nReads * cellMetrics.CROP_noGuide))) #Saturation 
        cellMetrics['CROP_assignedGuide'] = guideAssignment(cellGuideTab, thresh) #guide assignemt of any guides that are = max UMI value for that cell
    

    #optional filtering of cell metrics, guide UMI matrices and subsequent sample metrics quantitation
    if allcellsCsv is not None:
        #filter cellMetrics to only those in the RNA data that match the passing criteria in the RNA data
        rna = pd.read_csv(allcellsCsv)
        if 'i7_alias' in rna:
           rna = rna.rename(columns={"i7_alias": "i5_alias"}) #easy merging for downstream filtering in case a legacy i7 plate was used
        if 'I7_alias' in rna:
           rna = rna.rename(columns={"I7_alias": "i5_alias"}) 
        rna.index = pd.MultiIndex.from_frame(rna[['Ligation_alias', 'RT_alias', 'i5_alias']])
        rna_filtered = rna[rna['pass']]
        # I am only filtering on passing criteria here, without also filtering on i5_alias. The assumption is that for most workflows these should align
        cellMetrics_Filtered = rna_filtered.merge(cellMetrics, how='left', left_index=True, right_index=True)
        cellBarcodeConcat(references,cellMetrics_Filtered)
        cellMetrics_Filtered.fillna(0, inplace=True)
        cellMetrics_Filtered.drop(rna_filtered.columns, axis=1, inplace=True)
        cellMetrics_Filtered.loc[cellMetrics_Filtered['CROP_max'] < thresh, 'CROP_assignedGuide'] = 'NONE' #Set all non-assigned guide cells to 'NONE'
        
        cellMetrics_Filtered = cellMetrics_Filtered.reset_index()       
        cellMetrics_Filtered = cellMetrics_Filtered.set_index('Cell_Barcode')
        cellMetrics_Filtered.to_csv(outDir / "cellMetrics.filtered.csv")
        metrics.meanReadsPerCell = round(cellMetrics_Filtered['CROP_nReads'].mean(0),1)
        metrics.nUMIPerCell = round(cellMetrics_Filtered['CROP_nUMI'].median(0),1)
        metrics.nUMIPerCell = round(cellMetrics_Filtered['CROP_nUMI'].median(0),1)
        metrics.passingPercent = (cellMetrics_Filtered['CROP_max']>=thresh).mean()
        cellGuideTab_Filtered = rna_filtered.merge(cellGuideTab, how='left', left_index=True, right_index=True)
        cellBarcodeConcat(references,cellGuideTab_Filtered)
        cellGuideTab_Filtered.fillna(0, inplace=True)
        cellGuideTab_Filtered.drop(rna_filtered.columns, axis=1, inplace=True)
        cellGuideTab_noBC = cellGuideTab_Filtered.drop(['Cell_Barcode'], axis=1)
        metrics.meanUMIPerGuide = round((cellGuideTab_noBC.sum(0)).mean(),1)

        cellGuideTab_Filtered = cellGuideTab_Filtered.reset_index(drop=True)
        cellGuideTab_Filtered = cellGuideTab_Filtered.set_index('Cell_Barcode')
        #cellGuideTab_Filtered.drop(['Ligation_alias', 'RT_alias', 'i5_alias'], axis=1, inplace=True)
        cellGuideTab_Filtered.to_csv(outDir / "guideTab.filtered.csv")

        cellGuideReads_Filtered = rna_filtered.merge(cellGuideReads, how='left', left_index=True, right_index=True)
        cellBarcodeConcat(references,cellGuideReads_Filtered)
        cellGuideReads_Filtered.fillna(0, inplace=True)
        cellGuideReads_Filtered.drop(rna_filtered.columns, axis=1, inplace=True)
        cellGuideReads_Filtered = cellGuideReads_Filtered.reset_index(drop=True)
        cellGuideReads_Filtered = cellGuideReads_Filtered.set_index('Cell_Barcode')
        #cellGuideReads_Filtered.drop(['Ligation_alias', 'RT_alias', 'i5_alias'], axis=1, inplace=True)

        cellGuideReads_Filtered.to_csv(outDir / "guideReads.filtered.csv")
    else:
        cellMetrics.loc[cellMetrics['CROP_max'] < thresh, 'CROP_assignedGuide'] = 'NONE' #reassigns any guides who's max is below the threshold to 'NONE'
        metrics.meanReadsPerCell = round(cellMetrics['CROP_nReads'].mean(0),1)
        metrics.nUMIPerCell = round(cellMetrics['CROP_nUMI'].median(0),1)
        metrics.meanUMIPerGuide = round((cellGuideTab.sum(0)).mean(),1)
        metrics.passingPercent = (cellMetrics['CROP_max']>=thresh).mean()

    cellBarcodeConcat(references,cellGuideTab)
    cellBarcodeConcat(references,cellMetrics)
    cellBarcodeConcat(references,cellGuideReads)

    #handle cells with unassinged guides in the non-filtered matrix
    cellMetrics.fillna(0, inplace=True)
    cellMetrics.loc[cellMetrics['CROP_max'] < thresh, 'CROP_assignedGuide'] = 'NONE'
    cellMetrics = cellMetrics.reset_index()
    cellMetrics = cellMetrics.set_index('Cell_Barcode')
    cellMetrics.to_csv(outDir / "cellMetrics.csv")

    cellGuideTab = cellGuideTab.reset_index(drop=True)
    cellGuideTab = cellGuideTab.set_index('Cell_Barcode')
    #cellGuideTab.drop(['Ligation_alias', 'RT_alias', 'i5_alias'], axis=1, inplace=True)
    cellGuideTab.to_csv(outDir / "guideTab.csv")

    cellGuideReads = cellGuideReads.reset_index(drop=True)
    cellGuideReads = cellGuideReads.set_index('Cell_Barcode')
    #cellGuideReads.drop(['Ligation_alias', 'RT_alias', 'i5_alias'], axis=1, inplace=True)
    cellGuideReads.to_csv(outDir / "guideReads.csv")
    
    metrics.print(outDir / "metrics.csv")
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Count guide sequences matches from bcParser output')
    parser.add_argument('barcodesCsv', metavar='BARCODES.csv', type=Path,
                        help='bcParser per-sample barcode match output')
    parser.add_argument('--allcellsCsv', metavar='ALLCELLS.csv', type=Path,
                        help='optional nf-rna per-cell summary statistics output to filter based on pass status')
    parser.add_argument('--recountGuideTab', metavar='RECOUNTGUIDETAB.csv', type=Path,
                        help='optional unfiltered guidetab.csv to recount data with new RNA matrix')
    parser.add_argument('--recountCellMetrics', metavar='RECOUNTCELLMETRICS.csv', type=Path,
                        help='optional unfiltered guidetab.csv to recount data with new RNA matrix')
    parser.add_argument('--recountGuideReads', metavar='RECOUNTGUIDEREADS.csv', type=Path,
                        help='optional unfiltered guideReads.csv to recount data with new RNA matrix')
    parser.add_argument('--thresh', metavar='THRESH', nargs='?', const=3, type=int, default=3,
                        help='Theshold for the number of UMIs per guide in order to call a guide detected')
    parser.add_argument('--references', metavar='REFERENCES', type=Path,
                        help='Path to folder containing barcode whitelists')
    parser.add_argument('--outDir', type=Path,
                        help='Directory for outputs')
    args = parser.parse_args()
    args.outDir.mkdir(parents=True, exist_ok=True)
    main(args.barcodesCsv, args.allcellsCsv, args.recountGuideTab, args.recountCellMetrics, args.recountGuideReads, args.thresh, args.references, args.outDir)
