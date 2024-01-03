# Sample table

A sample table (e.g. [samples.csv](examples/samples.csv)) file is used to list the samples included in an analysis run, their sample barcode (RT) sequences and sample-specific analysis parameters.

It is a comma separated file (csv), with a header line (column names), followed by one sample per line. 
The first column is required to be `sample` and contains the name of each sample. These sample names need to be *UNIQUE*. This is especially important in the case of multiple final distribution plates, where the same sample as encoded by the original RT well will need to be found in multiple libraries and so their names need to be unique to maintain proper analysis downstream. The necessary column descriptions can be found below. _Column names are case sensitive!_

 Column | Description | Example
:---- | ---- | :----:
sample | Sample name | Foobar-2
barcodes | RT-plate wells used for this sample | 1A-2H
libName | Name for the overall sequencing library / fastq files | ScaleCROP
libIndex | Library barcode sequence corresponding to overall fastqName to separate out CROP reads| [BARCODE SEQUENCE]
guide | File name of guide sequence .tsv (each file should be in "guides" directory) | guides.txt
allCells | File name of ScaleRna output SAMPLE.allcells.csv corresponding the RNA libraries cell metrics (ultimately placed in allCells directory) | PBMC.allcells.csv

* `sample` and `libName` should consist only of letters, numbers, dash (`-`) and dot (`.`)
* A single `libName` should be used for an entire ScaleCROP sequencing library, not one `libName` per sample loaded into the RT plate.
    * In the case of multiple final distribution plates as in the Extended Throughput Kit, `libName` should correspond to each unique PLATE, and will have its own libIndex associated with each PLATE. See [samples.csv](examples/samples.csv) for illustration.
* When running from pre-existing fastq file input, `libName` should match the first part of the fastq file name for this sample, e.g.: `Foo1` for `Foo1_*.fastq.gz`.
* See [guides.txt](guides.md).
* See [allCells.csv](allCells.md).

## Demultiplexing samples
During analysis the sequencing data is first converted into library fastq files (`libName` column). If multiple samples were included in one sequencing library, these are then demultiplexed based on the sample (RT) barcodes. E.g.

sample | barcodes
-- | --
Foo | 1A-6H
Bar | 7A-12H

The RT wells used for each sample are given in `barcodes` as either
* An individual value (`1A`)
* A range of wells (`1A-2H`)
    * Wells are sorted first by number then by letter, i.e. `1A-1H`.
    * Note that all ranges are read in **column-wise** order; e.g. 1A-2C, refers to 1A-1H (all of column 1) plus 2A-2C.
* A list of values or ranges, separated by semicolon (`;`) (`1A;2A-2D`)
