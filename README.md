# ScaleBio Seq Suite: CRISPR Workflow

This is a Nextflow workflow to run analysis of ScaleBio single-cell CRISPR sequencing libraries. It processes data from sequencing reads to single-cell outputs (gene-count matrix, etc.), and QC reports.

## Getting started
* First install [Nextflow](http://www.nextflow.io) (22.04 or later)
* Download this workflow to your machine
* Install [dependencies](docs/dependencies.md)
* Create a [samples.csv](docs/samplesCsv.md) table for your samples
* Create [runParams.yml](docs/analysisParameters.md), specifying inputs and analysis options for your run
* Launch the workflow for your run

## Inputs
* Sequencing reads
    * Path to the Illumina Sequencer RunFolder (bcl files)
    * If you prefer to start from fastq files, generated outside (before) this workflow, see [Fastq generation](docs/fastqGeneration.md).
* samples.csv
    * A .csv file listing all samples in the analysis, optionally split by RT barcode. See [samples.csv](docs/samplesCsv.md).
* guides
    * Path to folder containing guide sequences to be detected and quantified by the workflow. See [guides.txt](docs/guides.md).
* allCells
    * Path to folder containing ScaleRna allCells.csv output files, one for each sample to be analyzed for the CROP analysis. See [allCells.csv](docs/allCells.md).
* outDir
    * Path to desired output directory for workflow analysis
## Outputs
The workflow produces per-sample QC reports (`metrics.csv`), a cell-by-guide umi count-matrix (`guideTab.filtered.csv`) and more; See [Outputs](docs/outputs.md) for a full list.


## Workflow Execution
### Workflow Example
Once the above inputs have been collected and dependancies installed, the workflow can be executed by the following:  

`nextflow run /PATH/TO/CROP -profile PROFILE -params-file /PATH/TO/runParams.yml --outDir output`

See [dependencies](docs/dependencies.md) for the best `PROFILE` to use on your system.

### Nextflow Command-line
**Note** that `nextflow` options are given with a single `-` (e.g. `-profile`), while workflow parameters (e.g. `--outDir`) are given with a double dash `--`.

See the [Nextflow command-line documentation](https://www.nextflow.io/docs/latest/cli.html) for the options to run `nextflow` on different systems (including HPC clusters and cloud compute).

## Configuration
### Specifying Analysis Parameters
Analysis parameters (inputs, options, etc.) can be defined either in a [runParams.yml](docs/examples/runParams.yml) file or directly on the nextflow command-line. See [analysisParameters](docs/analysisParameters.md) for details on the options.

### Config File
In addition to the analysis parameters, a user-specific nextflow configuration file can be used for system settings (compute and storage resources, resource limits, storage paths, etc.):

`-c path/to/user.config`

See [Nextflow configuration](https://www.nextflow.io/docs/latest/config.html) for the way different configuration files, parameter files and the command-line interact.

## Dependency Management
Different options to provide all required dependencies are described [here](docs/dependencies.md). Follow one approach there and then run nextflow with the corresponding `-profile`.

## Running in the cloud
Nextflow itself supports running using [AWS](https://www.nextflow.io/docs/latest/aws.html), [Azure](https://www.nextflow.io/docs/latest/azure.html) and [Google Cloud](https://www.nextflow.io/docs/latest/google.html). 

In addition [Nextflow tower](https://tower.nf) offers another simple way to manage and execute nextflow workflows in Amazon AWS.

# Versions and Updates
See the [Change log](changelog.md)

# License
For License information see [licenseNotice.md](docs/licenseNotice.md).

By purchasing product(s) and downloading the software product(s) of ScaleBio, You accept all of the terms of the [License Agreement](LICENSE.md). If You do not agree to these terms and conditions, You may not use or download any of the software product(s) of ScaleBio. 


