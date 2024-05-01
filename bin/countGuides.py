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
from scipy.sparse import csr_matrix, lil_matrix
from scipy.io import mmwrite


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

def create_sparse_matrix(data_dict, column_file):
    """Generate the sparse matrix for counts and UMI's from the dictionary output by countGuideReads"""

    # Read column names from TSV file
    with open(column_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        guides = [row[1] for row in reader]

    # Create mappings for row and column indices
    cell_mapping = {row: idx for idx, row in enumerate(data_dict.keys())}
    guide_mapping = {col: idx for idx, col in enumerate(guides)}

    # Initialize lists to store row indices, column indices, and data for the sparse matrix
    rows = []
    cols = []
    umi_data = []
    reads_data = []
    no_guide_values = {row_key: 0 for row_key in data_dict.keys()} 

    # Populate the lists with data from the dictionary
    for row_key in data_dict.keys():
        row_idx = cell_mapping[row_key]
        for guide in data_dict[row_key].guides:
            if guide in guide_mapping:
                col_idx = guide_mapping[guide]
                rows.append(row_idx)
                cols.append(col_idx)
                umi_data_val = data_dict[row_key].guides[guide].nUmis(1)
                umi_data.append(umi_data_val)

                reads_data_val = sum(data_dict[row_key].guides[guide].counts.values())
                reads_data.append(reads_data_val)
            elif guide == "noGuide":
                no_guide_values[row_key] = data_dict[row_key].guides[guide].nUmis(1)

    # Create the sparse matrix
    umi_sparse_matrix = csr_matrix((umi_data, (rows, cols)), shape=(len(cell_mapping), len(guide_mapping)))

    reads_sparse_matrix = csr_matrix((reads_data, (rows, cols)), shape=(len(cell_mapping), len(guide_mapping)))
    

     # Create a separate series for "noMatch" values mapped according to row_mapping
    no_guide_series = pd.Series({row_key: no_guide_values[row_key] for row_key in data_dict.keys()})

    no_guide_df = no_guide_series.to_frame(name='noGuide')

    return umi_sparse_matrix, reads_sparse_matrix, cell_mapping, guide_mapping, no_guide_df

def filter_sparse_matrix(sparse_matrix, MetaData, row_mapping):
    """Rebuild sparse matrix by only the cell barcodes that passed in the allCells.csv from RNA"""
    # Filter rows based on "pass" column in MetaData
    filtered_indices = MetaData.index[MetaData['pass']].tolist()

    # Initialize new row mapping and filtered sparse matrix
    filtered_row_mapping = {alias: idx for idx, alias in enumerate(filtered_indices)}
    filtered_sparse_matrix = lil_matrix((len(filtered_indices), sparse_matrix.shape[1]))

    # Populate filtered sparse matrix based on original row mapping and new filtered row mapping
    for alias, idx in filtered_row_mapping.items():
        if alias in row_mapping.keys():
            filtered_sparse_matrix[idx] = sparse_matrix[row_mapping[alias]]
    
    filtered_sparse_matrix = filtered_sparse_matrix.tocsr()

    return filtered_sparse_matrix, filtered_row_mapping

def find_second_highest_value_per_row(sparse_matrix):
    """Return value for the second highest detected guide per cell"""
    # Initialize an array to store second highest values
    second_highest_values = np.zeros(sparse_matrix.shape[0])

    # Iterate over rows
    for i in range(sparse_matrix.shape[0]):
        # Get row data
        row_data = sparse_matrix.data[sparse_matrix.indptr[i]:sparse_matrix.indptr[i + 1]]

        # Find the second highest value
        if len(row_data) > 1:
            row_data_sorted = np.sort(row_data)
            second_highest_values[i] = int(row_data_sorted[-2])

    return second_highest_values

def guide_assignment(sparse_matrix, threshold, column_mapping):
    """Creates a list of guides (columns) that passed the minimum UMI threshold within each cell"""
    # Initialize an empty list to store column names
    column_names_list = []

    # Iterate over rows
    for i in range(sparse_matrix.shape[0]):
        # Find column indices with values above threshold for the current row
        indices_above_threshold = np.where(sparse_matrix[i].A >= threshold)[1]

        # Map column indices to column keys using column mapping dictionary
        column_keys = [list(column_mapping.keys())[idx] for idx in indices_above_threshold]

        # Join column keys with ";" or append "Fail" if no columns pass
        if column_keys:
            column_names_list.append(";".join(column_keys))
        else:
            column_names_list.append("Fail")

    return column_names_list


def cellBarcodeConcat(referencePath,rows_dict):
    """Parse alieas in cellMetrics object per cell and create a series of the actual barcode sequences concatenated together like we do for RNA"""
    ligation_df = pd.read_csv(f'{referencePath}/3lvlRNA_lig.txt', sep="\t", names=['Barcode','Alias'])
    ligation_dict = dict(zip(ligation_df.Alias,ligation_df.Barcode))
    rt_df = pd.read_csv(f'{referencePath}/3lvlRNA_rt.txt', sep="\t", names=['Barcode','Alias'])
    rt_dict = dict(zip(rt_df.Alias,rt_df.Barcode))
    pcr_df = pd.read_csv(f'{referencePath}/3lvlRNA_pcr.txt', sep="\t", names=['Barcode','Alias'])
    pcr_dict = dict(zip(pcr_df.Alias,pcr_df.Barcode))
    alias_df = pd.DataFrame(rows_dict.keys(), columns = ['Ligation_alias', 'RT_alias', 'i5_alias'])
    alias_df['Ligation_alias'] = alias_df['Ligation_alias'].replace(ligation_dict)
    alias_df['RT_alias'] = alias_df['RT_alias'].replace(rt_dict)
    alias_df['i5_alias'] = alias_df['i5_alias'].replace(pcr_dict)
    alias_df['Cell_Barcode'] = alias_df['Ligation_alias'].astype(str) + alias_df['RT_alias'].astype(str) + alias_df['i5_alias'].astype(str)
    cell_barcode = alias_df['Cell_Barcode']
    return cell_barcode


def main(barcodesCsv: Path, allcellsCsv: Path, thresh: int, references: Path, guides_file: Path, outDir: Path):

    cellGuideReads, metrics, cellReads = countGuideReads(barcodesCsv)
    # Create  a table with UMI counts for each cell X guide combo
    # cellGuideCounts: Dict[]= {}
    cellMetrics = pd.DataFrame.from_dict(cellReads.counts, orient='index', columns=['CROP_nReads'])
    cellMetrics.index = pd.MultiIndex.from_tuples(cellMetrics.index, names=['Ligation_alias', 'RT_alias', 'i5_alias'])



    umi_sparse_counts_matrix, reads_sparse_counts_matrix, rows_mapped, cols_mapped, noGuides = create_sparse_matrix(cellGuideReads,guides_file)

    cell_barcodes = cellBarcodeConcat(referencePath=references, rows_dict=rows_mapped)

    cellMetrics['Cell_Barcode'] = cell_barcodes.values

    #Continue updating cellMetrics df

    cellMetrics['CROP_noGuide'] = noGuides['noGuide']/cellMetrics['CROP_nReads'] #Percentage of reads with no guide detected
    cellMetrics['CROP_nUMI'] = umi_sparse_counts_matrix.sum(axis=1) # number of unique guide UMI detected 
    cellMetrics['CROP_nGuide'] = umi_sparse_counts_matrix.getnnz(axis=1) # number of unique guide's detected, like unique genes detected
    cellMetrics['CROP_max'] = umi_sparse_counts_matrix.max(axis=1).toarray().flatten() # max number of UMI detected across the guides within a cell
    cellMetrics['CROP_second'] = find_second_highest_value_per_row(umi_sparse_counts_matrix) # second highest number of UMI of any guide within a cell
    cellMetrics['CROP_purity'] = (cellMetrics['CROP_max'] / cellMetrics['CROP_nUMI']) 
    cellMetrics['CROP_topTwo'] = ((cellMetrics['CROP_max']+cellMetrics['CROP_second']) / cellMetrics['CROP_nUMI'])
    cellMetrics['CROP_minorFrac'] = cellMetrics['CROP_second'] / cellMetrics['CROP_max']
    cellMetrics['CROP_guides'] = (umi_sparse_counts_matrix>=thresh).sum(1).A.flatten() # how many guides per cell passed a minimum detection threshold
    cellMetrics['CROP_Saturation'] = 1 - (cellMetrics.CROP_nUMI / (cellMetrics.CROP_nReads - (cellMetrics.CROP_nReads * cellMetrics.CROP_noGuide))) #Saturation 
    cellMetrics['CROP_assignedGuide'] = guide_assignment(umi_sparse_counts_matrix, thresh, cols_mapped) #guide assignemt of any guides that are = max UMI value for that cell
    
    cellMetrics = cellMetrics.reset_index()

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
        rna_filtered = rna_filtered.rename(columns = {'Unnamed: 0':'Cell_Barcode'})

        filtered_umi_sparse_matrix, filtered_row_mapping = filter_sparse_matrix(umi_sparse_counts_matrix, rna, rows_mapped)

        filtered_reads_sparse_matrix, _ = filter_sparse_matrix(reads_sparse_counts_matrix, rna, rows_mapped)


        filtered_cell_barcodes = cellBarcodeConcat(referencePath=references, rows_dict=filtered_row_mapping)

        # I am only filtering on passing criteria here, without also filtering on i5_alias. The assumption is that for most workflows these should align
        cellMetrics_Filtered = pd.merge(rna_filtered,cellMetrics, on='Cell_Barcode', how='left', suffixes=('_RNA','_ENRICH')).fillna(0)
        
        cellMetrics_Filtered.rename(columns=lambda x: x.replace('_ENRICH', ''), inplace=True)
        cellMetrics_Filtered = cellMetrics_Filtered[cellMetrics.columns]
        cellMetrics_Filtered.loc[cellMetrics_Filtered['CROP_max'] < thresh, 'CROP_assignedGuide'] = 'Fail' #Set all non-assigned guide cells to 'Fail'
        
        cellMetrics_Filtered = cellMetrics_Filtered.reset_index()       
        cellMetrics_Filtered = cellMetrics_Filtered.set_index('Cell_Barcode')
        cellMetrics_Filtered.to_csv(outDir / "cellMetrics.filtered.csv")

        
        metrics.meanReadsPerCell = round(cellMetrics_Filtered['CROP_nReads'].mean(0),1)
        metrics.nUMIPerCell = round(cellMetrics_Filtered['CROP_nUMI'].median(0),1)
        metrics.nUMIPerCell = round(cellMetrics_Filtered['CROP_nUMI'].median(0),1)
        metrics.passingPercent = (cellMetrics_Filtered['CROP_max']>=thresh).mean()
        metrics.meanUMIPerGuide = round((filtered_umi_sparse_matrix.sum(0)).mean(),1)

        #write out matrices
        filtered_umi_sparse_matrix = filtered_umi_sparse_matrix.transpose()
        mmwrite(f"{outDir}/filtered.umi.matrix.mtx", filtered_umi_sparse_matrix)
        filtered_cell_barcodes.to_csv(f"{outDir}/filtered.barcodes.tsv", sep='\t', index=False, header=False)

        filtered_reads_sparse_matrix = filtered_reads_sparse_matrix.transpose()
        mmwrite(f"{outDir}/filtered.reads.matrix.mtx", filtered_reads_sparse_matrix)

    else:
        cellMetrics.loc[cellMetrics['CROP_max'] < thresh, 'CROP_assignedGuide'] = 'NONE' #reassigns any guides who's max is below the threshold to 'NONE'
        metrics.meanReadsPerCell = round(cellMetrics['CROP_nReads'].mean(0),1)
        metrics.nUMIPerCell = round(cellMetrics['CROP_nUMI'].median(0),1)
        metrics.meanUMIPerGuide = round((umi_sparse_counts_matrix.sum(0)).mean(),1)
        metrics.passingPercent = (cellMetrics['CROP_max']>=thresh).mean()
    

    #handle cells with unassinged guides in the non-filtered matrix
    umi_sparse_counts_matrix = umi_sparse_counts_matrix.transpose()
    mmwrite(f"{outDir}/umi.matrix.mtx", umi_sparse_counts_matrix)
    cell_barcodes.to_csv(f"{outDir}/barcodes.tsv", sep='\t', index=False, header=False)

    reads_sparse_counts_matrix = reads_sparse_counts_matrix.transpose()
    mmwrite(f"{outDir}/reads.matrix.mtx", reads_sparse_counts_matrix)
    
    features_df = pd.DataFrame(cols_mapped.keys(), columns=['feature'])
    features_df.to_csv(f"{outDir}/features.tsv", sep='\t', index=False, header=False)

    cellMetrics = cellMetrics.set_index('Cell_Barcode')
    cellMetrics.to_csv(outDir / "cellMetrics.csv")

    metrics.print(outDir / "metrics.csv")
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Count guide sequences matches from bcParser output')
    parser.add_argument('barcodesCsv', metavar='BARCODES.csv', type=Path,
                        help='bcParser per-sample barcode match output')
    parser.add_argument('--allcellsCsv', metavar='ALLCELLS.csv', type=Path,
                        help='optional nf-rna per-cell summary statistics output to filter based on pass status')
    parser.add_argument('--thresh', metavar='THRESH', nargs='?', const=3, type=int, default=3,
                        help='Theshold for the number of UMIs per guide in order to call a guide detected')
    parser.add_argument('--references', metavar='REFERENCES', type=Path,
                        help='Path to folder containing barcode whitelists')
    parser.add_argument('--guides_file', metavar='GUIDES_FILE', type=Path,
                        help='Path to file containing guides list used for detection')
    parser.add_argument('--outDir', type=Path,
                        help='Directory for outputs')
    args = parser.parse_args()
    args.outDir.mkdir(parents=True, exist_ok=True)
    main(args.barcodesCsv, args.allcellsCsv, args.thresh, args.references, args.guides_file, args.outDir)
