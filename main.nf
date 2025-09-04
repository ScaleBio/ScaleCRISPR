nextflow.enable.dsl=2

include { INPUT_READS } from './modules/input_reads.nf'
include { shouldRunLibDetection } from './modules/utils.nf'
include { expandPath } from './modules/utils.nf'
include { reorderList; MERGING } from './modules/merging.nf'
include { SampleStats; CropAnalysis } from './modules/sample_report.nf'
include { resolveCrisprToRnaMapping } from './modules/utils.nf'

// Prepare samples.csv with defaults, rename legacy columns, etc.
process RegularizeSamplesCsv {
	label 'small'
	cache 'deep'

	// Publish user provided samples csv to output directory
	publishDir params.outDir, pattern: "samples.in.csv", mode: 'copy', saveAs: {'samples.csv'}

	input: 
	path("samples.in.csv")
	// Directory containing the library structure definition (barcode sequence lists are loaded from here)
	path("references")
	// Filename of the library structure definition .json
	val(libStructName)
	// Whether data was generated using QuantumScale assay
	val(quantum)
	path("guidesDir") // Directory containing guide sequences
	path(crisprToRnaMapping) // Crispr to RNA libName mapping file
	
	output: 
	path("samples.csv"), emit: samples_csv
	path("samples.in.csv")

	script:
	libStruct = "references/$libStructName"
	opts=""
	if (params.splitFastq) {
		opts += "--splitSample "
	}
	if (params.reporting) {
		opts += "--reporting "
	}
	if (params.resultDir) {
		opts += "--resultDir ${params.resultDir} "
	}
	if (quantum) {
		opts += "--quantum "
	}
	if (params.fastqDir) {
		opts += "--fastq "
	}
	if (crisprToRnaMapping) {
		opts += "--crisprToRnaMapping ${crisprToRnaMapping} "
	}
	"""
	regularize_samples_csv.py samples.in.csv --libraryStruct $libStruct --guide_dir guidesDir $opts
	"""
}

//Generate counts, reads, and UMI matrices
process CountGuides {
tag "$subsample"
label 'crispr'

input:
path(libStructDir)
val(libStructName)
val(file_ext)
tuple val(sample), val(subsample), path("ubam_*.${file_ext}"), path(guides) // seqTagFinder bam file with guide in tag and guides file

output:
tuple val(sample), val(subsample), path("${sample}.${subsample}.umi.matrix.mtx"), path("${sample}.${subsample}.reads.matrix.mtx"), path("${sample}.${subsample}.cellMetrics.csv"), emit: count_raw

script:
libStruct = "${libStructDir}/${libStructName}"
"""
	countGuides.py . --libStruct $libStruct --thresh ${params.thresh} --references ${libStructDir} --guides_file ${guides} --outDir ./ --id ${sample}.${subsample}
"""
}

process MergeCountRaw {
    tag "$sample"
	label 'large'

    input:
    tuple val(sample), path(umiMtx), path(readsMtx), path(cellMetrics)

    output:
    tuple val(sample), path("${sample}.umi.matrix.mtx"), path("${sample}.reads.matrix.mtx"), path("${sample}.cellMetrics.csv"), emit: countRaw

    script:
    """
    concat_columnar.py $cellMetrics --outputFile "${sample}.cellMetrics.csv"
    merge_mtx.py $umiMtx --outputFile "${sample}.umi.matrix.mtx"
    merge_mtx.py $readsMtx --outputFile "${sample}.reads.matrix.mtx"
    """
}

process FilterBarcodes {
tag "$sample"
publishDir { file(params.outDir) / "samples" / sample }, pattern: "filtered.reads.matrix.mtx", mode: 'copy', saveAs: { "filtered_reads_mtx/matrix.mtx" }
publishDir { file(params.outDir) / "samples" / sample }, pattern: "filtered.barcodes.tsv", mode: 'copy', saveAs: { "filtered_reads_mtx/barcodes.tsv" }
publishDir { file(params.outDir) / "samples" / sample }, pattern: "features.tsv", mode: 'copy', saveAs: { "filtered_reads_mtx/features.tsv" }
publishDir { file(params.outDir) / "samples" / sample }, pattern: "filtered.umi.matrix.mtx", mode: 'copy', saveAs: { "filtered_umi_mtx/matrix.mtx" }
publishDir { file(params.outDir) / "samples" / sample }, pattern: "filtered.barcodes.tsv", mode: 'copy', saveAs: { "filtered_umi_mtx/barcodes.tsv" }
publishDir { file(params.outDir) / "samples" / sample }, pattern: "features.tsv", mode: 'copy', saveAs: { "filtered_umi_mtx/features.tsv" }
publishDir { file(params.outDir) / "samples" / sample }, pattern: "reads.matrix.mtx", mode: 'copy', saveAs: { "raw_reads_mtx/matrix.mtx" }
publishDir { file(params.outDir) / "samples" / sample }, pattern: "barcodes.tsv", mode: 'copy', saveAs: { "raw_reads_mtx/barcodes.tsv" }
publishDir { file(params.outDir) / "samples" / sample }, pattern: "features.tsv", mode: 'copy', saveAs: { "raw_reads_mtx/features.tsv" }
publishDir { file(params.outDir) / "samples" / sample }, pattern: "umi.matrix.mtx", mode: 'copy', saveAs: { "raw_umi_mtx/matrix.mtx" }
publishDir { file(params.outDir) / "samples" / sample }, pattern: "barcodes.tsv", mode: 'copy', saveAs: { "raw_umi_mtx/barcodes.tsv" }
publishDir { file(params.outDir) / "samples" / sample }, pattern: "features.tsv", mode: 'copy', saveAs: { "raw_umi_mtx/features.tsv" }
publishDir { file(params.outDir) / "samples" / sample }, pattern: "cellMetrics.csv", mode: 'copy', saveAs: { "${sample}.cellMetrics.csv" }
publishDir { file(params.outDir) / "samples" / sample }, pattern: "cellMetrics.filtered.csv", mode: 'copy', saveAs: { "${sample}.cellMetrics.filtered.csv" }

input:
tuple val(sample), path("umi.matrix.mtx"), path("reads.matrix.mtx"), path("cellMetrics.csv"), path(barcodes), path("features.tsv")
path(crisprToRnaMapping) // Crispr to RNA libName mapping file

output:
tuple val(sample), path("filtered.umi.matrix.mtx"), path("cellMetrics.filtered.csv"), emit: countFiltered
tuple val(sample), path("filtered.reads.matrix.mtx"), emit: readsFiltered
tuple val(sample), path("umi.matrix.mtx", includeInputs: true)
tuple val(sample), path("reads.matrix.mtx", includeInputs: true)
tuple val(sample), path("cellMetrics.csv", includeInputs: true)
tuple val(sample), path("filtered.barcodes.tsv")
tuple val(sample), path("barcodes.tsv")
tuple val(sample), path("features.tsv", includeInputs: true)

script:
	if (crisprToRnaMapping) {
		opts = "--crisprToRnaMapping ${crisprToRnaMapping} "
	}
"""
	filter_barcodes.py --umiMatrix umi.matrix.mtx --readsMatrix reads.matrix.mtx \
		--thresh ${params.thresh} --barcodes $barcodes --cellMetrics cellMetrics.csv $opts
"""
}

process LibraryReport {
input:
	path(lib_template)
	tuple val(libs), path(metrics)
output:
	path("library.report.html")
publishDir { file(params.outDir) / "reports" / "library" }, mode: 'copy'
script:
"""
	generate_library_report.py --libs ${libs.join(" ")} --metrics $metrics --template $lib_template --outDir .
"""
}

workflow {
	ParamLogger.initialise(workflow, params, log)
	// Load rna library structure json
	libJson = expandPath(params.libStructure, file(projectDir) / "references")
	libJsonContents = Utils.loadJson(libJson)
	crisprToRnaMapping = resolveCrisprToRnaMapping(libJson)
	libraryInfo = [
		rnaLibraryStructureFile: libJson,
		rnaTrimAdapter: libJsonContents["trimAdapter"],
		rnaGenomicR1: libJsonContents.get('genomic_r1', false),
		rnaGenomicR2: libJsonContents.get('genomic_r2', false),
		quantum: libJsonContents.get('quantum', false),
		sampleBarcode: libJsonContents["sample_barcode"],
		crisprToRnaMappingFile: crisprToRnaMapping,
	]
	RegularizeSamplesCsv(
		file(params.samples, checkIfExists:true),
		libraryInfo.rnaLibraryStructureFile.getParent(),
		libraryInfo.rnaLibraryStructureFile.getName(),
		libraryInfo.quantum,
		file(params.guidesDir),
		libraryInfo.crisprToRnaMappingFile
	)
	samplesCsv = RegularizeSamplesCsv.out.samples_csv
	samples = samplesCsv.splitCsv(header:true, strip:true)
	samples.dump(tag:'samples')
	runLibDetection = shouldRunLibDetection(file(params.samples, checkIfExists:true), libraryInfo.quantum)
	// INPUT_READS runs bcParser on fastq files to demux per sample to unaligned bam files
	INPUT_READS(
		samples,
		samplesCsv,
		libraryInfo,
		params.runFolder,
		params.fastqDir,
		runLibDetection
	)

	file_ext = params.ultimaCramDir ? 'cram' : 'bam'
	INPUT_READS.out.ubam
		.groupTuple(by: [0,1])
		.combine(
			INPUT_READS.out.samples.map {
				[it.id, file(params.guidesDir) / it.guide]
			},
			by: 0
		).dump(tag: 'countGuidesInput')
		.set { countGuidesInput }

	CountGuides(libJson.getParent(), libJson.getName(), file_ext, countGuidesInput)
	CountGuides.out.count_raw
		.map { id, _subsample, rawUmi, readsUmi, cellMetrics ->
		// drop subsample that was used for splitting
		[id, rawUmi, readsUmi, cellMetrics]
		}
		.set { perSampleCountRaw }
	if (params.splitFastq) {
		perSampleCountRaw
			// collect all cellMetrics, barcodes, and umi.sparse.mtx for each lib and sample combination
			.groupTuple()
			| MergeCountRaw
			| set { perSampleCountRaw }
	}

	// resolve the result directory for each sample
	barcodes = INPUT_READS.out.samples.map {
		if (it.resultDir == null) {
			ParamLogger.throwError("Sample ${it.id} does not have a resultDir specified in samples.csv.")
		}
		barcodesTsv = file(it.resultDir).resolve("samples/${it.sample}.${it.rnaLibName}.filtered.matrix/barcodes.tsv.gz")
		if (!barcodesTsv.exists()) {
			barcodesTsv = file(it.resultDir).resolve("samples/${it.sample}_libraries/${it.sample}.${it.rnaLibName}.filtered.matrix/barcodes.tsv.gz")
		}
		if (!barcodesTsv.exists()) {
			ParamLogger.throwError("Sample ${it.sample}.${it.rnaLibName} was not found in the resultDir: ${it.resultDir}.")
		}
		[it.id, barcodesTsv]
	}

	FilterBarcodes(
		perSampleCountRaw
		    .join(barcodes)
			.join(
				INPUT_READS.out.samples.map {
					[it.id, file(params.guidesDir) / it.guide]
				}
			),
		libraryInfo.crisprToRnaMappingFile
	)

	SampleStats(
		perSampleCountRaw
			.map { sample, _rawUmi, _readsUmi, cellMetrics ->
				[sample, cellMetrics]
			}
			.join(FilterBarcodes.out.countFiltered
					.map { sample, _filteredUmi, filteredMetrics ->
						[sample, filteredMetrics]
					})
			.join(
				INPUT_READS.out.samples.map {
					[it.id, file(params.guidesDir) / it.guide]
				}
			)
	)

	LibraryReport(
		file(projectDir) / "templates" / "library_report_template.html",
		// collect all libraries and their metrics, e.g. [['QSR-1', 'QSR-2'], [/path/to/QSR-1.metrics.json, /path/to/QSR-2.metrics.json]]
		INPUT_READS.out.metrics
			.collect(flat: false)
			.transpose()
			.collect(flat: false)
			.map {
				libs, metrics ->
				[
					reorderList(libs, libs),
					reorderList(metrics, libs)
				]
			}
	)

	CropAnalysis(
		libJson.name,
		file(projectDir) / "templates" / "qc_report_template.html",
		FilterBarcodes.out.countFiltered
		 	.join(SampleStats.out.sample_metrics)
			.join(
				INPUT_READS.out.samples.map { [it.id, it.barcodes] }
			)
	)

	if (params.merge) {
		MERGING(
			INPUT_READS.out.samples,
			libJson,
			FilterBarcodes.out.countFiltered
				.map { sample, _filteredUmi, filteredMetrics ->
					[sample, filteredMetrics]
				},
			perSampleCountRaw
				.map { sample, _rawUmi, _readsUmi, cellMetrics ->
					[sample, cellMetrics]
				},
			FilterBarcodes.out.readsFiltered,
			FilterBarcodes.out.countFiltered
				.map { sample, filteredUmi, _filteredMetrics ->
					[sample, filteredUmi]
				},
			INPUT_READS.out.samples.map {
				[it.id, file(params.guidesDir) / it.guide]
			},
			barcodes,
		)
	}

}

//// When the workflow completes, publish 'workflow_info.json' containing information on pipeline metadata
workflow.onComplete {
	testing = false
	def workflow_data = ["Nextflow Information":
		["Execution status": "${ workflow.success ? 'OK' : 'failed' }",
		"Run name": "$workflow.runName",
		"Pipeline completion timestamp": "$workflow.complete",
		"Git repo URL": "$workflow.repository",
		"Configuration files": "$workflow.configFiles",
		"Container": "$workflow.container",
		"Command line executed": "$workflow.commandLine",
		"Configuration profile": "$workflow.profile",
		"Start timestamp": "$workflow.start",
		"Stop timestamp": "$workflow.complete",
		"Exit status": "$workflow.exitStatus",
		"Error message": "$workflow.errorMessage",
		"Revision": "$workflow.revision",
		"Commit ID": "$workflow.commitId",
		"Duration": "$workflow.duration"]
	]
	def params_data = ["Parameters": [:]]
	for (p in params) {
		if (!p.key.contains('-')) {
			params_data."Parameters".put("$p.key", "$p.value")
		}
    }
	def manifest_data = ["Workflow Manifest": [:]]
	for (p in workflow.manifest.getProperties()) {
		p = p.toString()
		def split_str = p.split("=")
		if (split_str[0].equals("name") || split_str[0].equals("version")) {
			manifest_data["Workflow Manifest"].put(split_str[0], split_str[1])
		}
	}

	def json_str = groovy.json.JsonOutput.toJson(manifest_data + workflow_data + params_data)
	def json_beauty = groovy.json.JsonOutput.prettyPrint(json_str)
	def workflow_info = file(params.outputDir) / "workflow_info.json"
	workflow_info.write(json_beauty)
}
