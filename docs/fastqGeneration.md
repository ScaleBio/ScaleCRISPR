# Fastq generation
The workflow can be started from a sequencer runFolder (bcl) [`--runFolder`]. In that case fastq files will be generated internally, using Illumina `bcl-convert`.

Alternatively it is also possible to generate fastq files upstream, for example when the ScaleCROP library is multiplexed together with other libraries during sequencing at a core facility or with a ScaleRNA library. In this case we recommend using [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html).

An example [samplesheet.csv](examples/samplesheet.csv) with typical options is included. Here all 96 i5 barcode sequences from the PCR plate are merged into one set of fastq files. If an index1 (i7) read is used to demultiplex the ScaleBio CRISPR library with other libraries in the sequencing run, usage of the `index` column should contain the constant i7 sequence of the CRISPR library, which can be found here in forward orientation: [3lvlCRISPR_p7.txt](../references/3lvlCRISPR_p7.txt). This becomes especially relevant when working with multiple final distribution plates. In this case, each sample name for each row in the samples.csv needs to corresopnd to each sample found PER PLATE uniquely.

## Index reads
For ScaleBio RNA and CRISPR libraries the RT and ligation barcodes are included in read 1, the PCR plate barcode is in the index1 read, while the PCR well barcode is in the index2 read. Hence we need to tell `bcl-convert` to generate index read fastqs using the `samplesheet.csv` setting: \
`CreateFastqForIndexReads,1`

# Using pre-generated fastq files as workflow input
Set `--fastqDir` to the directory containing the fastq files. 
The file names should follow the pattern `<LibraryName>_..._<Read>_...fastq.gz`, where
* `Name` is the library name ( set in the `libName` column in `samples.csv`)
* `Read` is one of `R1`, `R2`, `I1`, `I2`
