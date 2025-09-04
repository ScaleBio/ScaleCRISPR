#!/usr/bin/env python
"""
Pre-process samples.csv to simplify parsing downstream.
Rename deprecriated column names and fill in defaults
Basic validation

Writes the normalized samples.csv to stdout
"""
import argparse
import csv
import sys
import pandas as pd
from pathlib import Path
from typing import Optional, Dict
from scale_utils.validation import validateName
from scale_utils.lib_json_parser import LibJsonParser


COL_NAMES_WHITELIST = [
    "id",
    "sample",
    "barcodes",
    "libName",
    "libIndex",
    "libIndex2",
    "split",
    "group",
    "resultDir",
    "rnaId",
    "expectedCells",
    "guide",
    "rnaLibName",
]


def makeUniqueSampleIds(cols, rows, quantum):
    """
    Add an 'id' column with a unique value for each row

    The same sampleName can occur more than once per library. (e.g. multiple PCR plates).
    Hence we can use "sampleName.libName" as unique id
    """
    for r in rows:
        validateName(r[cols.index("sample")], display_name="sample")
        validateName(r[cols.index("libName")], display_name="libName", other_chars="-")
        if quantum and "libIndex2" in cols:
            r.insert(0, f"{r[cols.index('sample')]}.{r[cols.index('libIndex2')]}")
        else:
            r.insert(0, f"{r[cols.index('sample')]}.{r[cols.index('libName')]}")
    # add column name at end so indices are correct in line above
    cols.insert(0, "id")


def fail(msg: str, sep=" "):
    """
    Helper function to print an error message and exit

    Args:
        msg: Error message
        sep: Separator to use when printing the error message
    """
    print(msg, file=sys.stderr, sep=sep)
    sys.exit(1)


def unwrapSemiColonSeparatedValues(cols: list, rows: list, index: int) -> list:
    """
    Unwrap each entry separated by a semi-colon into its own row

    Args:
        cols: Column headers
        rows: List of rows from the csv
        index: Index of the column to unwrap

    Returns:
        final_rows: List of rows with unwrapped values
    """
    final_rows = []
    for r in rows:
        # Unwrap each entry separated by a semi-colon into its own row
        if ";" in r[index]:
            split_libs = r[index].split(";")
            split_libs = [libname.strip() for libname in split_libs]
            for idx, split_libname in enumerate(split_libs):
                unwrapped_row = r[:]
                unwrapped_row[index] = split_libname
                final_rows.append(unwrapped_row)
        else:
            final_rows.append(r)
    return final_rows


def readSamplesCsv(samplesCsv: Path, quantum: bool, fastq_input: bool, reporting: bool) -> tuple[list, list[list[str]]]:
    """
    Load samples.csv into a list.

    Perform basic validation and cleanup.
    Args:
        samplesCsv: Path to samples.csv file
        quantum: Flag to indicate analysis is for data generated on a ScaleQuantum kit
        fastq_input: Flag to indicate that fastq files will be used as input to the workflow
        reporting: Flag to indicate this is a reporting run

    Returns:
        cols: Column headers
        rows: List of rows from the csv
    """
    rows = []
    # Create a mapping of lowercase column names to their original names
    col_names_mapping = {col.lower(): col for col in COL_NAMES_WHITELIST}
    # utf8-sig ignores 'BOM' characters at the beginning of the file, which some users might
    # accidentally add and which will break downstream parsing of samples.csv column headers.
    with open(samplesCsv, encoding="utf-8-sig") as csvfile:
        samples = csv.reader(csvfile)
        cols = next(samples)  # CSV column header (in same order as values in rows)
        # Normalize column names
        cols = [
            col_names_mapping[col.lower()] if col.lower() in col_names_mapping else fail(f"Invalid column name {col}")
            for col in cols
        ]
        # Trim empty columns
        while not cols[-1]:
            cols = cols[:-1]
        if len(cols) != len(set(cols)):
            fail("Duplicate column names, check samplesheet")
        if "sample" not in cols:
            fail("'sample' column missing (the sample name)")
        for row in samples:
            if not row or len(row[0].strip()) == 0:
                continue
            # Remove leading and trailing spaces from all values in the row
            row = [row_element.strip() for row_element in row]
            # Trim empty columns
            while not row[-1] and len(row) > len(cols):
                row = row[:-1]

            if len(row) != len(cols):
                fail(f"Unexpected number of columns: {row}", sep="\n")
            rows.append(row)
    if "libName" not in cols and "libIndex" in cols:
        for r in rows:
            if ";" in r[cols.index("libIndex")]:
                raise ValueError("libName is required if libIndex lists multiple index sequences")
    if "libName" in cols and "libIndex" in cols:
        for r in rows:
            if ";" in r[cols.index("libIndex")] and ";" in r[cols.index("libName")]:
                raise ValueError("Cannot list multiple libIndex if listing multiple libName")
    if quantum and not reporting:
        if "libIndex2" in cols and "libName" in cols:
            raise ValueError("Please provide a value for only the libIndex2 column or the libName column, not both")
        if not fastq_input:
            if "libName" in cols:
                raise ValueError("When starting from bcl files, do not provide a value for the libName column")
    # Normalize index sequences
    for idx_col in ["libIndex", "libIndex2"]:
        if idx_col in cols:
            for r in rows:
                index = "".join(r[cols.index(idx_col)].split())  # Remove whitespace
                index = index.upper()
                index = ";".join(sorted(index.split(";")))  # sort sequences
                r[cols.index(idx_col)] = index
    final_rows = []
    if quantum and "libName" in cols:
        final_rows = unwrapSemiColonSeparatedValues(cols, rows, cols.index("libName"))
    else:
        final_rows = rows
    if quantum and "libIndex2" in cols:
        final_rows = unwrapSemiColonSeparatedValues(cols, final_rows, cols.index("libIndex2"))
    else:
        final_rows = final_rows
    if "guide" not in cols:
        raise ValueError("Please provide a value for the guide column")
    return cols, final_rows


def check_whether_barcode_ranges_overlap(all_samples_barcode_range: Dict, sample_barcode_fname: Path):
    """
    Check whether the user supplied barcode ranges overlap amongst samples
    Throw exception if there is an overlap

    Args:
        all_samples_barcode_range: Dictionary of barcode ranges for all samples with libName as key
        sample_barcode_fname: Path to sample barcode whitelist file
    """
    barcode_whitelist = pd.read_csv(sample_barcode_fname, sep="\t", names=["well", "barcode"], index_col=False)
    for libName in all_samples_barcode_range:
        # Needs to be rest for each libName
        verbose_barcodes_list = []
        for sample_barcodes in all_samples_barcode_range[libName]:
            # String can be a single barcode or a range of barcodes separated by semi-colon
            for semi_colon_separated in sample_barcodes.split(";"):
                if "-" in semi_colon_separated:
                    # Get well coordinate that corresponds to starting of the barcodes range for a sample
                    starting = barcode_whitelist.index[barcode_whitelist["well"] == semi_colon_separated.split("-")[0]][
                        0
                    ]
                    # Get well coordinate that corresponds to end of the barcodes range for a sample
                    end = barcode_whitelist.index[barcode_whitelist["well"] == semi_colon_separated.split("-")[1]][0]
                    # Retrieve all the well coordinates that correspond to the barcodes range for a sample
                    all_barcodes = barcode_whitelist.loc[starting:end, "well"].tolist()
                    # extend because all_barcodes is a list
                    verbose_barcodes_list.extend(all_barcodes)
                else:
                    verbose_barcodes_list.append(semi_colon_separated)
        # Check whether the barcode ranges overlap amongst individual samples
        if len(verbose_barcodes_list) != len(set(verbose_barcodes_list)):
            fail("The barcodes range mentioned for each sample overlap amongst individual samples")


def write_unwrapped_barcodes_range_samples_csv(cols: list, rows: list, sample_barcode_fname: Path):
    """
    Write a csv with a singular well coordinate for the barcodes column on each line with same id,
    sample name and library name
    Relevant for analyzing ultima data where each well coordinate maps to a single file

    Args:
        cols: Column headers
        rows: Rows from the csv
        sample_barcode_fname: Path to sample barcode whitelist file
    """
    unwrapped_rows = []
    barcode_whitelist_wells = pd.read_csv(sample_barcode_fname, sep="\t", names=["well", "barcode"])["well"].to_list()
    for r in rows:
        unwrapped_rows.append(
            unwrap_barcodes_range(
                r[cols.index("barcodes")], r[cols.index("id")], r[cols.index("libName")], barcode_whitelist_wells
            )
        )
    flattened_unwrapped_rows = [
        item for sublist in unwrapped_rows for item in (sublist if isinstance(sublist, list) else [sublist])
    ]
    with open("barcode_range_unwrapped_samples.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["id", "libName", "barcode"])
        writer.writerows(flattened_unwrapped_rows)


def get_all_barcodes_in_range(barcode_range: str, barcode_whitelist_wells):
    """
    Get all barcodes in the range defined by @barcode_range

    Args:
        barcode_range: Range of barcodes
        barcode_whitelist_wells: All barcode well coordinates in whitelist file

    Returns:
        List of all barcodes in the range
    """
    starting_barcode = barcode_range.split("-")[0]
    ending_barcode = barcode_range.split("-")[1]
    starting_index = barcode_whitelist_wells.index(starting_barcode)
    ending_index = barcode_whitelist_wells.index(ending_barcode)
    return barcode_whitelist_wells[starting_index : ending_index + 1]


def unwrap_barcodes_range(
    barcodes_range: str, sample_id: str, libName: str, barcode_whitelist_wells: list
) -> list[tuple[str, str]]:
    """
    Unwrap the barcodes range into individual barcodes

    Args:
        barcodes_range: Range of barcodes
        sample_id: Sample id
        libName: Library name
        barcode_whitelist_wells: All barcode well coordinates in whitelist file

    Returns:
        List of tuples containing sample_id, libName, barcode
    """
    if ";" in barcodes_range:
        unwrapped_barcodes = []
        split_barcodes = barcodes_range.split(";")
        for split_barcode in split_barcodes:
            if "-" in split_barcode:
                unwrapped_barcodes.extend(get_all_barcodes_in_range(split_barcode, barcode_whitelist_wells))
            else:
                unwrapped_barcodes.append(split_barcode)
    elif "-" in barcodes_range:
        unwrapped_barcodes = get_all_barcodes_in_range(barcodes_range, barcode_whitelist_wells)
    else:
        return (sample_id, libName, barcodes_range)
    return [(sample_id, libName, barcode) for barcode in unwrapped_barcodes]


def add_all_alias_for_index(
    cols: list[str],
    rows: list[list[str]],
    all_aliases: list[str],
    col_name: str,
) -> list[list[str]]:
    """
    For a particular barcode (corresponds to @read), add all aliases from the whitelist file for each sample

    Args:
        cols: Column headers
        rows: Rows from the csv
        all_aliases: All aliases from whitelist file for the given read
        col_name: Column name to add aliases for

    Returns:
        List of rows with all aliases added
    """
    rows_with_all_aliases = []

    for r in rows:
        for alias in all_aliases:
            r_copy = r[:]
            r_copy.insert(cols.index(col_name), alias)
            rows_with_all_aliases.append(r_copy)

    return rows_with_all_aliases


def check_guide_files_exist(cols: list, rows: list, guide_dir: Path):
    """
    Check whether the guide files exist in the guide directory

    Args:
        cols: Column headers
        rows: Rows from the csv
        guide_dir: Path to the directory containing the guide whitelist files
    """
    for r in rows:
        guide_file = r[cols.index("guide")]
        if not (guide_dir / guide_file).exists():
            fail(f"Guide file {guide_file} does not exist in {guide_dir}")


def main(
    samplesCsv: Path,
    splitFastq: bool,
    reporting: bool,
    libraryStruct: Path,
    quantum: bool,
    resultDir: Optional[str],
    samples_csv_fname: str,
    fastq_input: str,
    guide_dir: Path,
    crispr_to_rna_mapping: Optional[dict[str, str]],
):
    """
    Writes normalized samples.csv to stdout
    Args:
        samplesCsv: Path to samples.csv file
        splitFastq: Split fastq files per RT well
        reporting: Generate a samples.csv for a merging/reporting only workflow
        quantum: Analysis is for data generated on a ScaleQuantum kit
        resultDir: Path or URL to previous workflow outputs
        libraryStruct: Path to library structure json
        samples_csv_fname: Name of the samples.csv file
        fastq_input: Flag to indicate that fastq files will be used as input to the workflow
        guide_dir: Path to the directory containing the guide whitelist files
        crispr_to_rna_mapping: Optional mapping from CRISPR to RNA libNames
    """
    lib_json_obj = LibJsonParser(libraryStruct)
    sample_barcode_whitelist_first_entry, sample_barcode_whitelist_second_entry = lib_json_obj.get_barcodes_range_of_whitelist_file(lib_json_obj.sample_barcode_fname)
    barcodes_range = f"{sample_barcode_whitelist_first_entry}-{sample_barcode_whitelist_second_entry}"
    # Load column headers and all lines from the csv into lists
    cols, rows = readSamplesCsv(samplesCsv, quantum, fastq_input, reporting)
    check_guide_files_exist(cols, rows, guide_dir)
    # Process csv (add defaults, rename columns etc.)
    # Default libName based on index
    if "libName" not in cols:
        if quantum:
            if "libIndex2" not in cols:
                cols.insert(len(cols), "libIndex2")
                rows = add_all_alias_for_index(
                    cols,
                    rows,
                    lib_json_obj.all_whitelist_contents_by_read["index2"],
                    "libIndex2",
                )
        if crispr_to_rna_mapping and "rnaLibName" not in cols and "libIndex2" in cols:
            cols.append("rnaLibName")
            for r in rows:
                r.append(crispr_to_rna_mapping.get(r[cols.index("libIndex2")]))
        # Add libName column after libIndex or at end of cols
        idx = cols.index("libIndex") if "libIndex" in cols else len(cols)
        for r in rows:
            if "libIndex" in cols:
                name = r[cols.index("libIndex")]
            elif "libIndex2" in cols and quantum:
                name = r[cols.index("libIndex2")]
            else:
                name = "ScaleRNA"
            r.insert(idx, name)
        # insert column name at idx after row updates
        cols.insert(idx, "libName")

    # Resolve sample names and IDs
    if cols[0] != "id":
        if "id" in cols:
            fail("id has to be the first column")
        else:
            makeUniqueSampleIds(cols, rows, quantum)
    ids = set()
    for r in rows:
        id = r[0]
        if id in ids:
            fail(f"Duplicate sample name: {id}")
        else:
            ids.add(id)

    all_samples_barcode_range = {}

    if not reporting:
        libNameIndex = cols.index("libName")
        for r in rows:
            all_samples_barcode_range[r[libNameIndex]] = []
        if "barcodes" not in cols:
            cols.append("barcodes")
            for r in rows:
                r.append(barcodes_range)
                all_samples_barcode_range[r[libNameIndex]].append(barcodes_range)
        else:  # check whether each entry for barcodes is empty or not
            barcodesIndex = cols.index("barcodes")
            for r in rows:
                if r[barcodesIndex].strip() == "":
                    r[barcodesIndex] = barcodes_range
                    all_samples_barcode_range[r[libNameIndex]].append(barcodes_range)
                else:
                    r[barcodesIndex] = r[barcodesIndex].replace(" ", "")
                    all_samples_barcode_range[r[libNameIndex]].append(r[barcodesIndex])
        check_whether_barcode_ranges_overlap(all_samples_barcode_range, lib_json_obj.sample_barcode_fname)
        if splitFastq:
            if "split" not in cols:
                cols.append("split")
                for r in rows:
                    r.append("true")
    # "group" is used for merging samples with --merge; default is by "sample"
    # (e.g. different libraries of the same sample)
    if "group" not in cols:
        nameInd = cols.index("sample")
        cols.append("group")
        for r in rows:
            r.append(r[nameInd])
    if resultDir and "resultDir" not in cols:  # Set a single resultDir for all samples
        cols.append("resultDir")
        for r in rows:
            r.append(resultDir)
    with open(samples_csv_fname, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        w.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Standardize samples.csv for the workflow (defaults, column names, etc."
    )
    parser.add_argument(
        "samples", metavar="SAMPLES.csv", help="CSV with sample names and information for the workflow run"
    )
    parser.add_argument(
        "--splitSample",
        help="Flag that denotes whether bcParser will split per RT well to ensure parallelization",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--reporting", help="set for use in merging/reporting only workflow", action="store_true", default=False
    )
    parser.add_argument(
        "--resultDir", help="Previous pipeline output directory (used as default for 'resultDir column)"
    )
    parser.add_argument("--libraryStruct", help="Library structure json file", type=Path)
    parser.add_argument(
        "--quantum",
        help="Flag to indicate analysis is for data generated on a ScaleQuantum kit",
        action="store_true",
        default=False,
    )
    parser.add_argument("--samples_csv_fname", help="Name of the samples.csv file", type=str, default="samples.csv")
    parser.add_argument(
        "--fastq",
        help="Flag to indicate that fastq files will be used as input to the workflow",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--guide_dir",
        help="Path to the directory containing the guide whitelist files",
        type=Path,
        required=True,
    )
    parser.add_argument("--crisprToRnaMapping", type=Path, help="Path to crispr to rna mapping file")
    args = parser.parse_args()
    
    crispr_to_rna_mapping = {}
    if args.crisprToRnaMapping:
        df = pd.read_csv(args.crisprToRnaMapping, sep="\t", header=None)
        crispr_to_rna_mapping = dict(zip(df[0], df[1]))
    
    main(
        samplesCsv=args.samples,
        splitFastq=args.splitSample,
        reporting=args.reporting,
        libraryStruct=args.libraryStruct,
        quantum=args.quantum,
        resultDir=args.resultDir,
        samples_csv_fname=args.samples_csv_fname,
        fastq_input=args.fastq,
        guide_dir=args.guide_dir,
        crispr_to_rna_mapping=crispr_to_rna_mapping,
    )