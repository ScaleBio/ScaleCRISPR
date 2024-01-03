# Analysis parameters

All analysis parameters can be set in a [runParams.yml](examples/runParams.yml), which is then passed to nextflow with `nextflow run -params-file runParams.yml`. 
Alternatively each option is this file can also be set on the nextflow command-line directly, overwriting the value in the parameter file. E.g.
`nextflow run --samples=samples.foo.csv`

*Note* that `nextflow` options such as `resume` or `dump-channels` are given with a single `-`, while workflow parameters (e.g. `samples`) are given with a double dash `--`.


## Inputs
The workflow can either start from an Illumina sequencer runFolder (bcl files) or a directory with fastq files. Specify either
* runFolder : path/to/runFolder <br>
OR
* fastqDir : path/to/fastqs

where fastqDir is a directory containing all input fastq files. See [Fastq Generation](fastqGeneration.md) for details on file names, etc.

When starting from a sequencer run folder the workflow uses Illumina [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) for automatic fastq generation.

### Sample Information
* samples : samples.csv

A [file](examples/samples.csv) listing all samples in the analysis with their names, barcode sequences and optional sample settings

### ScaleRna allCells
* allCellsDir : path/to/allCells

Path to folder containing ScaleRna allCells.csv output files, one for each sample to be analyzed for the CROP analysis. See [allCells.csv](allCells.md). This is used to ensure that the CRISPR output matrix has the same features as the RNA matrix. Each file is referenced inthe samples.csv. An easy way to to use this is to point to the "samples" output directory of the ScaleRNA workflow, which will contain the necessary files for the paired RNA portion of the library.

### Guide Reference Sheets
* guidesDir : path/to/guides

Path to folder containing guide sequences to be detected and quantified by the workflow. See [guides.txt](guides.md). These are in tab-separated format (.tsv) and the constituent file names are referenced in [samples.csv](examples/samples.csv)
    * IMPORTANT NOTE: These sequences must all be the same lenght for the pipeline to proceed, and if they are not the software will trim them for you down to the lowest sequence length



## Optional and Advanced parameters
### Threshold
* thresh : minThresh

Integer value corresponding to minimum threshold (>= thresh) of guide UMI detection within a cell for that guide to be passing and thus "assigned" to that cell. Defaults to 3 within the workflow, but can be adjusted based on your data.

Run `nextflow run path/to/ScaleCRISPR --help` for a description of available options and see the example [runParams.yml](examples/runParams.yml) file.

System options (compute resource requirements, etc.) as well as all parameter defaults, are in the workflow [nextflow.config](../nextflow.config).

### Library Structure Definition
* libStructure : libV1.json
* libTemplate : libGuideSeq.TEMPLATE.json

The libStructure JSON file defines 
* Where in the reads cell-barcodes are found
* What the list of allowed barcode sequences is

The default file, for our standard product configuration, is included in [references/libV1.json](../references/libV1.json) and a per library file will be generated for each sample based on a template for guide detection within the workflow [references/libGuideSeq.TEMPLATE.json](../references/libGuideSeq.TEMPLATE.json). These two files should match in the exception of the "guide" barcode portion in the template file, which has placeholder values that will be calculated for and repalced within the workflow for proper handling per sample / guide reference sheet. 
