#!/usr/bin/env python

import csv
import gzip
import json
from Bio import SeqIO
import argparse
from dataclasses import dataclass
import numpy as np
import pandas as pd
from pathlib import Path

@dataclass
class guideMetrics:
    """Overall read statistics for the dataset"""
    guideLength: int = 0 # guide sequence length 
    startPosition: int = 0 # start position of the guide sequence on the read

    def print(self, outFn: Path):
        with open(outFn, 'w') as out:
            print("Length", f"{self.guideLength:}", sep=',', file=out)
            print("StartPosition", f"{self.startPosition:}", sep=',', file=out)



def qc_guide_input(guide_list, outDir):
    # Read the TSV file
    outputMetrics = guideMetrics()
    with open(guide_list, 'r', newline='') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        rows = list(reader)

    # Check if the lengths are the same
    lengths = set(len(row[0]) for row in rows)
    min_length = min(lengths)
    outputMetrics.guideLength = min_length
    if len(lengths) > 1:
        
        #print(f"The minimum length of the strings is {min_length}")

        # Trim strings to the minimum length
        trimmed_rows = [[row[0][:min_length], row[1]] for row in rows]
        print(f"Trimmed strings to the minimum length `{min_length}`")
        with open(outDir / "guides.output.txt", 'w', newline='') as output_file:
            writer = csv.writer(output_file, delimiter='\t')
            writer.writerows(trimmed_rows)
        finalGuide = (outDir / "guides.output.txt")
    else:
        print(f"All guides are the same length: {min_length}, no trimming necessary")
        trimmed_rows = rows
        with open(outDir / "guides.output.txt", 'w', newline='') as output_file:
            writer = csv.writer(output_file, delimiter='\t')
            writer.writerows(trimmed_rows)
        finalGuide = (outDir / "guides.output.txt")
    # Check for duplicate string values
    string_count = {}
    duplicate_strings = []
    for row in trimmed_rows:
        print(row)
        string = row[0]
        name = row[1]
        if string in string_count:
            if string not in duplicate_strings:
                duplicate_strings.append(string)
        else:
            string_count[string] = name
    #write out duplicate guides to file
    if len(duplicate_strings) > 0:
        print("Duplicate string values found:")
        with open(outDir / "Duplicate.guides.txt", "w") as file:
            for string in duplicate_strings:
                line = f"Sequence: {string}, Name: {string_count[string]}"
                file.write(line)
    else:
        print("No duplicate string values found.")
    return trimmed_rows,outputMetrics, finalGuide

def calculate_average_position(fastq_file, target_sequences, maxreads):
    results = {}
    linesRead = 0
    with gzip.open(fastq_file,"rt") as handle: # unzip fastq file
        for record in SeqIO.parse(handle, "fastq"): 
            sequence = str(record.seq) #create string of read sequence
            linesRead += 1
            if linesRead == maxreads:
                break #Break when maxreads have been reached
            for target_sequence, _ in target_sequences: #read in guide sequences
                if target_sequence in sequence: #check match of guide sequence to read sequence
                    start_position = sequence.index(target_sequence) # record position of match if present
                    #keep track of guide sequence position, or create new entry if first
                    if target_sequence in results: 
                        results[target_sequence].append(start_position)
                    else:
                        results[target_sequence] = [start_position]

        averages = {}
        for target_sequence, positions in results.items():
            average_position = sum(positions) / len(positions)
            averages[target_sequence] = average_position

        return averages

#Recursive function for value replacement in dictionary. To be used to update json file template to new values
def dict_replace_value(d, old, new):
    x = {}
    for k, v in d.items():
        if isinstance(v, dict):
            v = dict_replace_value(v, old, new)
        elif isinstance(v, list):
            v = list_replace_value(v, old, new)
        elif isinstance(v, str):
            if v == old:
                if isinstance(new, int):
                    v = int(new)
                else:
                    v = new
        x[k] = v
    return x

#Recursive function for value replacement in dictionary. To be used to update json file template to new values
def list_replace_value(l, old, new):
    x = []
    for e in l:
        if isinstance(e, list):
            e = list_replace_value(e, old, new)
        elif isinstance(e, dict):
            e = dict_replace_value(e, old, new)
        elif isinstance(e, str):
            if e == old:
                if isinstance(new, int):
                    e = int(new)
                else:
                    e = new
        x.append(e)
    return x



# reads in a template libGuideSeq.json file and replaces placeholder values with supplied new values, whereever it finds matches
def update_libGuideSeq_json(libGuideSeq,start,length, outDir, guidesPath):
    with open(libGuideSeq, 'r') as file:
        template = json.load(file)
    updated = dict_replace_value(template,'START_REPLACE', start)
    updated = dict_replace_value(updated,'LENGTH_REPLACE', length)
    updated = dict_replace_value(updated,'PATH_REPLACE', "guides.output.txt")
    with open(outDir / 'updated_LibGuideSeq.json', 'w') as file:
        json.dump(updated, file, indent=4)



def main(fastq: Path, guides: Path, libGuideTemplate: Path, maxreads: int, outDir: Path):
    guides_qc, outputMetrics, finalGuidePath = qc_guide_input(guides, outDir)
    target_sequences = []
    for row in guides_qc:
        target_sequence = row[0]
        target_name = row[1]
        target_sequences.append((target_sequence, target_name))

    average_positions = calculate_average_position(fastq, target_sequences, maxreads)
    overall_average = round(np.mean(list(average_positions.values())))
    outputMetrics.startPosition = overall_average
    outputMetrics.print(outDir / "guideQCMetrics.csv")
    with open(outDir / "averagePosition.txt", "w") as file:
        header = f"Average position of all sequences: {overall_average}\n"
        file.write(header)
        for target_sequence, average_position in average_positions.items():
            line = f"Average position of '{target_sequence}': {average_position}\n"
            file.write(line)
    if libGuideTemplate is not None: # Optional libGuideSeq.json update from template based on calculated position and length
        update_libGuideSeq_json(libGuideTemplate,outputMetrics.startPosition,outputMetrics.guideLength, outDir, finalGuidePath)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculate average start BP postion of guide sequences in R2 of sequencing data for bcParser input')
    parser.add_argument('--fastq', metavar='FASTQ.gz', type=Path,
                        help='R2 sequencing data of a CROP-seq library')
    parser.add_argument('--guides', metavar='GUIDES.txt', type=Path,
                        help='tsv file of guide sequences that will ultimately be queried by bcParser')
    parser.add_argument('--libGuideTemplate', metavar='LIBGUIDETEMPLATE.json', type=Path,
                        help='template file of libGuideSeq.json to be updated for proper bcParser processing downstream')
    parser.add_argument('--maxreads', metavar='MAXREADS', nargs='?', const=100000, type=int, default=100000,
                        help='The number of reads the script will search through for average calculation')
    parser.add_argument('--outDir', type=Path,
                        help='Directory for outputs')
    args = parser.parse_args()
    args.outDir.mkdir(parents=True, exist_ok=True)
    main(args.fastq, args.guides, args.libGuideTemplate, args.maxreads, args.outDir)