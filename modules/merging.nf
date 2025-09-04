include { SampleStats; CropAnalysis } from './sample_report.nf'

// Reorder list based on what would sort refList
// In practice, refList would be sample names and the other lists would be workdir paths
def reorderList(list, refList) {
    def argsort = refList.sort(false).collect { name -> refList.indexOf(name) } // false means don't sort in-place
    argsort.collect { index -> list[index] }
}

// Merge cell metrics and filtered mtx for individual libraries
// Assumes different cells in each subsample
process MergeMetrics {
	tag "$group"
    publishDir { file(params.outDir) / "samples" / group }, pattern: "reads.mtx", mode: 'copy', saveAs: { "filtered_reads_mtx/matrix.mtx" }
    publishDir { file(params.outDir) / "samples" / group }, pattern: "filtered.barcodes.tsv", mode: 'copy', saveAs: { "filtered_reads_mtx/barcodes.tsv" }
    publishDir { file(params.outDir) / "samples" / group }, pattern: "features.tsv", mode: 'copy', saveAs: { "filtered_reads_mtx/features.tsv" }
    publishDir { file(params.outDir) / "samples" / group }, pattern: "umi.mtx", mode: 'copy', saveAs: { "filtered_umi_mtx/matrix.mtx" }
    publishDir { file(params.outDir) / "samples" / group }, pattern: "filtered.barcodes.tsv", mode: 'copy', saveAs: { "filtered_umi_mtx/barcodes.tsv" }
    publishDir { file(params.outDir) / "samples" / group }, pattern: "features.tsv", mode: 'copy', saveAs: { "filtered_umi_mtx/features.tsv" }
    publishDir { file(params.outDir) / "samples" / group }, pattern: "cellMetrics.filtered.csv", mode: 'copy', saveAs: { "${group}.cellMetrics.filtered.csv" }

	input:
	tuple val(group), path("cellMetrics_???.csv"), path("features.tsv"), path("reads_???.mtx"), path("umi_???.mtx"), path("barcodes_???.tsv.gz")

    output:
    tuple val(group), path("umi.mtx"), path("cellMetrics.filtered.csv"), emit: merged
    tuple val(group), path("reads.mtx")
    tuple val(group), path("features.tsv", includeInputs: true)
    tuple val(group), path("filtered.barcodes.tsv")
	
	script:
	"""
    concat_columnar.py cellMetrics_*.csv --outputFile cellMetrics.filtered.csv
    concat_columnar.py barcodes_*.tsv.gz --outputFile filtered.barcodes.tsv
    concat_mtx.py reads_*.mtx --id $group --threads $task.cpus --memory '$task.memory' --outputFile reads.mtx
    concat_mtx.py umi_*.mtx --id $group --threads $task.cpus --memory '$task.memory' --outputFile umi.mtx
	"""
}

workflow MERGING {
take:
    samples // each row of samples.csv parsed into a channel 
    libJson // library structure json file
    cellMetricsFiltered // Cell metrics for passing cells
    cellMetrics // Metrics for all cell-barcode combinations
    filteredReadsMtx // Read counts matrix for guides
    filteredUmiMtx // Umi counts matrix for guides
    guides // Path to guides files
    barcodes // barcodes.tsv.gz from RNA passing cells
main:
    samples
        .map { tuple(it.id, it.group) }
        .join(cellMetricsFiltered)
        .join(guides)
        .join(filteredReadsMtx)
        .join(filteredUmiMtx)
        .join(barcodes)
        .groupTuple(by:1) // All samples grouped for merging by library/group
        .filter{ it[0].size() > 1 } // Only merge samples with multiple libraries
        .dump(tag:'mergeSamplesInput')
        .map {
            ids, group, passingCellMetrics, guidesFiles, readsMtx, umiMtx, barcodesTsv ->
            // Compute indices that would sort the samples list for deterministic output
            [
                group + ".merged",
                reorderList(passingCellMetrics, ids),
                guidesFiles.first(),
                reorderList(readsMtx, ids),
                reorderList(umiMtx, ids),
                reorderList(barcodesTsv, ids)
            ]
        }
        .dump(tag:'mergeMetricsInput')
        .set{ mergeMetricsInput }
    MergeMetrics(mergeMetricsInput)

    samples
        .map { tuple(it.id, it.group) }
        .join(cellMetrics)
        .join(cellMetricsFiltered)
        .join(guides)
        .groupTuple(by:1) // All samples grouped for merging by library/group
        .filter{ it[0].size() > 1 } // Only merge samples with multiple libraries
        .map {
            ids, group, rawCellMetrics, passingCellMetrics, guidesFiles ->
            [
                group + ".merged",
                reorderList(rawCellMetrics, ids),
                reorderList(passingCellMetrics, ids),
                guidesFiles.first() // Implicit assumption that guides are the same for all samples in a group
            ]
        }
        .dump(tag:'sampleStatsInput')
        .set{ sampleStatsInput }
        SampleStats(sampleStatsInput)

    samples
        .map { tuple(it.id, it.group, it.barcodes)}
        .groupTuple(by:1) // All samples grouped for merging by library/group
        .filter{ it[0].size() > 1 } // Only merge samples with multiple libraries
        .map {
            _ids, group, barcodesRange ->
            [
                group + ".merged",
                Utils.combineSampleBarcodes(barcodesRange)
            ]
        }
        .dump(tag:'mergeBarcodesInput')
        .set { mergedBarcodesRange }
    MergeMetrics.out.merged
        .join(SampleStats.out.sample_metrics)
        .map {
            group, umiMtx, passingCellMetrics, metrics ->
            [
                group,
                umiMtx,
                passingCellMetrics,
                metrics,
            ]
        }
        .join(mergedBarcodesRange)
        .dump(tag:'mergeAnalysisInput')
        .set{ mergeAnalysisInput }

        CropAnalysis(
            libJson.name,
            file(projectDir) / "templates" / "qc_report_template.html",
            mergeAnalysisInput,
        )
}