#!/usr/bin/env python


import argparse
import pysam
import json
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Tuple
from collections import defaultdict
import csv
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, lil_matrix
from scipy.io import mmwrite


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


def countGuideReads(bam_dir: Path=Path(".")):
    """Count reads for each cell X guide X UMI combo from bcParser bam output"""
    # BC -> (GUIDE,UMI) -> ReadCount
    # BC is a tuple of three strings for the three cell-barcode levels, e.g. (RT, Ligation, PCR) for 3lvl
    cellGuideReads: Dict[Tuple[str, ...], GuideReadCounts] = defaultdict(GuideReadCounts)
    cellReads = cellReadCounts()
    bam_files = list(bam_dir.glob("*.bam")) + list(bam_dir.glob("*.cram"))
    for bam in bam_files:
        with pysam.AlignmentFile(bam, "rb", check_sq=False) as samfile:
            for read in samfile.fetch(until_eof=True):
                umi = read.get_tag("UM")
                cell = tuple(read.get_tag("CB").split("+"))
                cellReads.addCount(cell)  # Count a read for the cell, having had a correct match to the barcode info
                if read.has_tag("tg"):
                    cellGuideReads[cell].addRead(read.get_tag("tg"), umi)
                else:
                    # Count reads that do not have a guide sequence but did have a correct cellular barcode
                    cellGuideReads[cell].addRead("noGuide", read.query_name)

    return cellGuideReads, cellReads

def create_sparse_matrix(counts_data_dict, column_file):
    """
    Generate the sparse matrix for counts and UMI's from the dictionary output by countGuideReads BC -> (GUIDE,UMI) -> ReadCount
    Args:
        counts_data_dict: dictionary output by countGuideReads in the format of BC -> (GUIDE,UMI) -> ReadCount
        column_file: the guide whitelist used by bcParser for guide detection, first column is barcode sequences and second column are the aliases
    Returns:
        umi_sparse_matrix: csr sparse matrix of umi counts, in cells x guides orientation
        reads_sparse_matrix: csr sparse matrix of read counts, in cells x guides orientation
        cell_mapping: dictionary of mapping of cell barocdes to their index with tuple of cell barcode alias's as keys, and corresponding index in matrix as values (rows here)
        guide_mapping: dictionary of mapping of guide alias to their index in the matrix (columns here)
        no_guide_df: single column pd.Dataframe containing the number of reads with no guide per cell. Eventually appends to cellMetrics and converted to a percentage
    """

    # Read in guide whitelist used for detection in bcParser to parse into a list that will become our columns in the matrix being built
    with open(column_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        guides = [row[0] for row in reader]

    # Create mappings for row and column indices
    cell_mapping = {row: idx for idx, row in enumerate(counts_data_dict.keys())}
    guide_mapping = {col: idx for idx, col in enumerate(guides)}

    # Initialize lists to store row indices, column indices, and data for the sparse matrix
    rows = []
    cols = []
    umi_data = []
    reads_data = []
    no_guide_values = {row_key: 0 for row_key in counts_data_dict.keys()} 

    # Populate the lists with data from the dictionary generated by countGuideReads, BC -> (GUIDE,UMI) -> ReadCount
    for row_key in counts_data_dict.keys():
        row_idx = cell_mapping[row_key]
        for guide in counts_data_dict[row_key].guides:
            if guide in guide_mapping:
                col_idx = guide_mapping[guide]
                rows.append(row_idx)
                cols.append(col_idx)
                umi_data_val = counts_data_dict[row_key].guides[guide].nUmis(1)
                umi_data.append(umi_data_val)

                reads_data_val = sum(counts_data_dict[row_key].guides[guide].counts.values())
                reads_data.append(reads_data_val)
            elif guide == "noGuide":
                no_guide_values[row_key] = counts_data_dict[row_key].guides[guide].nUmis(1)

    # Create the sparse matrix, cells x guides
    umi_sparse_matrix = csr_matrix((umi_data, (rows, cols)), shape=(len(cell_mapping), len(guide_mapping)))

    reads_sparse_matrix = csr_matrix((reads_data, (rows, cols)), shape=(len(cell_mapping), len(guide_mapping)))
    

    # Create a separate series for "noMatch" values mapped according to row_mapping
    no_guide_series = pd.Series({row_key: no_guide_values[row_key] for row_key in counts_data_dict.keys()})

    no_guide_df = no_guide_series.to_frame(name='noGuide')

    return umi_sparse_matrix, reads_sparse_matrix, cell_mapping, guide_mapping, no_guide_df


def find_second_highest_value_per_row(sparse_matrix):
    """Return value for the second highest detected guide per cell"""
    # Initialize an array to store second highest values
    second_highest_values = np.zeros(sparse_matrix.shape[0])

    # Iterate over rows
    for i in range(sparse_matrix.shape[0]):
        # Get row data
        row_data = sparse_matrix.getrow(i).data

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
        indices_above_threshold = (sparse_matrix[i, :] >= threshold).nonzero()[1]

        # Map column indices to column keys using column mapping dictionary
        column_keys = [list(column_mapping.keys())[idx] for idx in indices_above_threshold]

        # Join column keys with ";" or append "Fail" if no columns pass
        if column_keys:
            column_names_list.append(";".join(column_keys))
        else:
            column_names_list.append("Fail")

    return column_names_list


def read_lib_json(lib_json: Path):
    with open(lib_json, 'r') as f:
        lib_struct = json.load(f)
    barcodes = []
    for bc in lib_struct['barcodes']:
        if bc.get("name") not in ["umi", "guide"]:
            barcodes.append(bc.get("alias"))
    return barcodes


def main(bam_dir: Path, lib_struct: Path, thresh: int, references: Path, guides_file: Path, outDir: Path, id: str):

    barcodes = read_lib_json(lib_struct) 
    # 1. Convert barcodesCsv to bam (which will be input, including multiple, in new version)
    # 2. Load library json to get mapping file, cell barcodes, and guide whitelists
    cellGuideReads, cellReads = countGuideReads(bam_dir)
    # Create  a table with UMI counts for each cell X guide combo
    cellMetrics = pd.DataFrame.from_dict(cellReads.counts, orient='index', columns=['CROP_nReads'])
    cellMetrics.index = pd.MultiIndex.from_tuples(cellMetrics.index, names=barcodes)

    umi_sparse_counts_matrix, reads_sparse_counts_matrix, _rows_mapped, cols_mapped, noGuides = create_sparse_matrix(cellGuideReads,guides_file)

    cellMetrics["Cell_Barcode"] = cellMetrics.index.map("+".join)

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

    #transpose to be guides x cells to be consistent with rna expression data
    umi_sparse_counts_matrix = umi_sparse_counts_matrix.transpose()
    mmwrite(f"{outDir}/{id}.umi.matrix.mtx", umi_sparse_counts_matrix)

    reads_sparse_counts_matrix = reads_sparse_counts_matrix.transpose()
    mmwrite(f"{outDir}/{id}.reads.matrix.mtx", reads_sparse_counts_matrix)
    
    cellMetrics = cellMetrics.set_index('Cell_Barcode')
    cellMetrics.to_csv(outDir / f"{id}.cellMetrics.csv")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Count guide sequences matches from bcParser output')
    parser.add_argument('bamDir', default='.', metavar='/path/ubam_1.bam', type=Path,
                        help='Path to seqTagFinder output bam(s) with guide alias in tag')
    parser.add_argument("--libStruct", type=Path, help="library structure json used for bcParser")
    parser.add_argument('--thresh', metavar='THRESH', nargs='?', const=3, type=int, default=3,
                        help='Theshold for the number of UMIs per guide in order to call a guide detected')
    parser.add_argument('--references', metavar='REFERENCES', type=Path,
                        help='Path to folder containing barcode whitelists')
    parser.add_argument('--guides_file', metavar='GUIDES_FILE', type=Path,
                        help='Path to file containing guides list used for detection')
    parser.add_argument('--outDir', type=Path,
                        help='Directory for outputs')
    parser.add_argument('--id', type=str, help='Sample id for output files')
    args = parser.parse_args()
    args.outDir.mkdir(parents=True, exist_ok=True)
    main(args.bamDir, args.libStruct, args.thresh, args.references, args.guides_file, args.outDir, args.id)
