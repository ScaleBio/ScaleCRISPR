# Outputs

Analysis workflow results for all samples in the run will appear in the output directory defined by `--outDir` (`ScaleCROP.out` by default). 


## Key output files
| Directory | File | Description |
|-----------|------|-------------|
| `reports` | `multiqc_report.html` | [MultiQC](https://multiqc.info/) report for fastq generation, fastQC and trimming |
| `samples` | `<sample>` | Per-sample outputs and analysis of the CRISPR workflow, described in detail below |
| `fastq` | `fastqc/*_fastqc.html` | [fastqc](https://github.com/s-andrews/FastQC) report for each fastq file in the sequencing library |
| | `Reports/` | Fastq generation summary reports from [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) |
| `barcodes/<libName>.demux` | `<sample>.fastq.gz` | Sample fastq files (Demultiplexed and barcode error-corrected); only included with `--fastqOut true` |


Within the output directory, there will then be per-sample (corresponding to the sample column in samples.csv) in the `samples` sub-directory that contains the sample specific output results. Within each, there is a a `<sample>.sample.metrics.csv` summary file, a `csv` directory containing output tables, a `references` folder containing the library strucutre reference and post-qc guides list used for detection, and an `analysis` folder containing figures output from the workflow. The majority of figures and analyses present correspond to the "filtered" tables, that contain only those cell barcodes that "passed" in the ScaleRna workflow. In the crispr workflow, cells that passed in RNA but did not have crispr library data will be padded with zeros. Unfiltered tables are also output, but these only correspond to cell barcodes found within the crispr library and have not been padded. The full results can be best summarized by the following:

* guideTab.filtered.csv : cells x guide UMI count matrix, with rows as cell index combos and columns as guides labeled as taken from the guides.txt

* guideReads.filtered.csv : cells x guide reads count matrix, with rows as cell index combos and columns as guides labeled as taken from the guides.txt

* SAMPLE.sample.metrics.csv : Per sample QC metrics regarding library performance

| Measure | ​Description​ |
|---------|-------------|
| Reads | Total number of reads |
| noGuide​ | % of reads that have no guide​ |
| Correct​ | % of usable reads |​
| ReadsPerCell | ​Mean reads per cell​ |
| nUMIPerCell​ | Median UMI per cell​ |
| meanUMIPerGuide | mean UMI per guide across all cells total |
| passingPercent | % of cells that had a passing Guide detected​ |

* cellMetrics.filtered.csv : Per cell QC metrics regarding UMI guide detection. The majority of these columns have the word CROP appended as a prefix so as not to cause conflit with potential merging with the allCells.csv for inclusion of metadata in downstream analysis.

| Measure | ​Description​ |
|---------|-------------|
| Ligation_alias | ligation plate well designation |
| RT_alias | RT plate well designation |
| PCR_alias | PCR or Final Distribution Plate well Designation |
| noGuide | % of reads with no guide detected​ |
| nUMI​ | Count of UMI detected​ |
| nGuide​ | Number of unique guides detected​ |
| Max​ | Top guide UMI count​ |
| Second​ | Second Highest guide UMI count​ |
| Purity​ | Proportion of UMI coming from max guide |​
| TopTwo​ | Proportion of UMI coming from both max and second |​
| MinorFrac​ | Ratio of second to max​ |
| Guides​ | Number of guides >= minium detection threshold​ |
| Saturation​ | Ratio of UMI detected over Usable Reads |
| assignedGuide​ | Guide assignment label, for any guides that passed the threshold​ |


There are a host of figures that are generated in addition to the output tables, and visualize the cellMetrics output calculated from the UMI counts matrix. These can be found in the `analysis` subdirectory of each sample, and a brief description of each output can be found below:

* Cell_Count_Containging_Passing_Guides.png : Histogram of cell numbers containing X number of "passing" guides as determined by the minimum theshold (default 3 UMI per guide). Zero corresponds to the number of cells that had no passing guides, One corresponds to the number of cells that had exactly one passing guide, so on and so forth for 2,3,4, etc...

* Guide_Purity_Density.png : Distribution of the Purity measures for each cell in the analysis

* Guide_vs_RNA_Counts_Scatter.png : Scatterplot of Top Guide UMI's vs Total RNA UMI's

* Max_Guide_UMI_Density.png : Distribution of guide UMI counts of the highest detected guide per cell

* nUMI_Distribution.png : Distribution of total guide UMI counts detected per cell

* nUMI_vs_nReads_Scatter.png : Scatterplot of Guide UMI's vs CRISPR library reads per cell

* Passing_Guide_and_Cell_Percent.png : The X axis encodes potential thresholds for minimum guide UMI detection PER guide that can be used to determine whether a guide is "passing" or not. The y axis encodes the percentage of features (cells or guides) that can pass. The Blue line corresponds to visualizing guide dropout, as we increase the minimum detection threshold, what percentage of guides have at least one cell in which that guide "passes" the threshold. The orange line corresonds to what percentage of cells have at least one guide passing at that threshold, and the green line describes the percentage of cells that have at least two guides passing the threshold. Finally, the dotted red line is placed at the threshold that was used for the analysis for determination of passing guides (default 3)

* Passng_Cell_Tresholds.txt : A tab delimited text file containing the passing percentage data plotted in the Passing_Guide_and_Cell_Percent.png figure

* Passing_Guide_Cell_Histogram.png : Histogram of the count of guides that passed in X number of cells

* Proportion_of_Max_Guide.png : ECDF plot representing the proportion of cells falling below the max guide value along the x axis

* Saturation_vs_nReads_Scatter.png : Scatter plot of guide UMI saturation vs total CRISPR library reads per cell, labeled by whether that cell had a guide assigned to that cell or not

* Guide_Cell_Numbers.txt : a tab delimited text file containing for each guide the number of cells that passed and were assigned. In the case of a cell being assigned multiple guides, each guide is counted independently. This can translate to cells being double counted, as we are summarizing in how many cells each guide was assigned.