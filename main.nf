nextflow.enable.dsl=2

import groovy.json.JsonOutput
import groovy.json.JsonSlurper

include { inputReads } from './modules/inputReads'

// Load .json file into object
// (reference genomes, library design, etc.)
def loadJson(json) {
	def jsonSlurper  = new JsonSlurper()
	return jsonSlurper.parse( json )
}

// Create a 'file()' from string 'path'
// 'path' could be null, a s3 URL, an absolute file-path or relative to 'baseDir'
def expandPath(path, baseDir) {
	if (path == null) { return null}
	if (path =~ /^s3/) { return file(path)}
	return baseDir.resolve(path)
}

// parse sample sheet into sample name and path to allcells from that run

process regularizeSamplesCsv {
input: 
	path("samples.in.csv")
output: 
	path("samples.csv")
publishDir "${params.outDir}", mode: 'copy'
cache 'deep'
label 'small'
script:
	opts=""
	if (params.splitFastq) {
		if (params.runFolder == null) {
			opts += "--splitSample"
		}
	}
"""
	regularizeSamplesCsv.py samples.in.csv $opts > samples.csv
"""
}


// this needs to operate on only the R2 of the reads folder
process qcGuideinput {
input:
	path(guidesDir)
	tuple(val(sample),path(fastq),path(barcode),val(guides))
	path(libStructDir)
	path(libGuideTemplate)
output:
	tuple(val(sample),path(fastq),path(barcode),path("tmpReferences/updated_LibGuideSeq.json"),path("tmpReferences/guides.output.txt"))
publishDir "${params.outDir}/samples/${sample}/references/", pattern: "*/updated_LibGuideSeq.json", saveAs:{"updated_LibGuideSeq.json"}, mode: 'copy'
publishDir "${params.outDir}/samples/${sample}/references/", pattern: "*/guides.output.txt", saveAs:{"post.qc.guides.txt"}, mode: 'copy'
label 'crispr'
tag "$sample"
script:
"""
	qcGuideInput.py --fastq ${fastq} --guides ${guidesDir}/${guides} --outDir ./tmpReferences --libGuideTemplate ${libStructDir}/${libGuideTemplate}
"""
}

// bcParser to find guide sequences
process barcode_assignment {
input:
	path(libStructDir)
	tuple(val(sample),path(fastq),path(barcode),path("tmpReferences/updated_LibGuideSeq.json"),path("tmpReferences/guides.output.txt"))
output:
	tuple(val(sample),path("${sample}/${sample}_*.barcodes.csv"), emit: barcodesCsv)
	tuple(val(sample),path("${sample}/metrics.json"), emit: bc_metrics)
tag "$sample"
script:
"""
	cp $libStructDir/* ./tmpReferences/
	bc_parser --lib-struct tmpReferences/updated_LibGuideSeq.json -v --read1 ${barcode} --read2 ${fastq} --out ${sample} --barcode-info --sample ${sample}
"""
}

//Generate counts, reads, and UMI matrices
process count_guides {
input:
	path(allCellsDir)
	path(libStructDir)
	tuple(val(sample), path(barcodesCsv),val(allcells)) // bcparser output barcodes.csv and corresponding sample allCells
output:
	tuple(val(sample), path("guideTab.filtered.csv"), path("cellMetrics.filtered.csv"),val(allcells), emit: count_filtered)
	tuple(val(sample), path("guideTab.csv"), path("cellMetrics.csv"),val(allcells), emit: count_raw)
	tuple(val(sample), path("guideReads.csv"), path("guideReads.filtered.csv"),val(allcells), emit: reads)
	tuple(val(sample),path("metrics.csv"), emit: count_metrics)
publishDir ("${params.outDir}/samples/${sample}/", pattern: "*.csv", mode: 'copy', saveAs: {fn -> 
	if (fn.endsWith("metrics.csv")) { "${sample}.sample.metrics.csv" }
	else { "csv/${fn}" }
	}
)

//publishDir "${params.outDir}/samples/${sample}/", pattern: "metrics.csv", saveAs:{"${sample}.sample.metrics.csv"}, mode: 'copy'
//publishDir "${params.outDir}/samples/${sample}/csv/", pattern: "*.csv", mode: 'copy'
label 'crispr'
tag "$sample"
script:
"""
	countGuides.py ${barcodesCsv} --allcellsCsv ${allCellsDir}/${allcells} --thresh ${params.thresh} --references ${libStructDir} --outDir ./
"""
}

//figure generation
process crop_analysis {
input:
	path(allCellsDir)
	tuple(val(sample), path("guideTab"), path("cellMetrics"),val("allcells")) // count_guides filtered cell x guide umi matrix and cell guide metrics, nf-rna output allcells.csv path
output:
	path("*.png")
	path("*.txt")// for future report generation
publishDir "${params.outDir}/samples/${sample}/analysis/", mode: 'copy'
label 'crispr'
tag "$sample"
script:
"""
	cropAnalysis.py --sample ${sample} --guideTabCsv ${guideTab} --allcellsCsv ${allCellsDir}/${allcells} --guideMetricsCsv ${cellMetrics} --thresh ${params.thresh} --outDir ./
"""
}

workflow {
	// be sure to join outputs on sample to ensure that files are properly paired
	// Call initialise function in Helper.groovy to pretty print important
	// parameters along with some additional information about the pipeline

	regularizeSamplesCsv(file(params.samples))
	samplesCsv = regularizeSamplesCsv.out
	samples = samplesCsv.splitCsv(header:true, strip:true)

	//3lvl libJson without the inclusion of the CROP guide detection, for the initial bcParser demux by sample performed by inputReads below
	libJson = expandPath(params.libStructure, file("${projectDir}/references/"))
	//call input reads and create per-sample fastq files split by PCR well. 
	inputReads(samples, samplesCsv, libJson, params.runFolder, params.fastqDir, params.fastqSamplesheet, params.trimFastq, params.fastqc, params.outDir, params.bclConvertParams, params.splitFastq, params.trimAdapt)

	fqsSorted = inputReads.out.fqs.toSortedList().flatMap()
	fqsSorted.dump(tag:'fqsSorted')
	alignFqJobGroups = fqsSorted.groupTuple()

	// Sort transcript and barcode files
	alignFqJobGroups = alignFqJobGroups.map {
		tuple(it[0], it[1].sort(), it[2].sort())
	}

	//generate input channels for QC of guide sequences
	guides = samples.map{
		tuple(it.sample, it.guide)
	}

	fastq_qc_ch = alignFqJobGroups.join(guides)
	fastq_qc_ch.dump(tag:'fastq_qc_ch')

	libTemplate = expandPath(params.libTemplate, file("${projectDir}/references/"))

	//qc guide input and write out proper libstruct.json for each sample
	qcGuideinput(params.guidesDir,fastq_qc_ch,libJson.getParent(),libTemplate)
	

	//perform bcParser detection of guide sequences with updated_libstructJson as output from qcGuide_input and fastq files from inputreads
	barcode_assignment(libJson.getParent(),qcGuideinput.out)

	//generate input channel of 
	allcells = samples.map{[it.sample, it.allCells]}

	count_input_ch = barcode_assignment.out.barcodesCsv.join(allcells)
	count_guides(params.allCellsDir,libJson.getParent(),count_input_ch)

	crop_analysis(params.allCellsDir,count_guides.out.count_filtered)

}

// Function that gets invoked when workflow completes
// Publishes a file named workflow_info.json that contains information about the pipeline execution
workflow.onComplete {
	testing = false
	def data = ["Workflow Information": ["Execution status": "${ workflow.success ? 'OK' : 'failed' }",
										 "Pipeline completion timestamp": "$workflow.complete",
										 "Git repo URL": "$workflow.repository",
										 "Configuration files": "$workflow.configFiles",
		 		 						 "Container": "$workflow.container",
				 						 "Command line executed": "$workflow.commandLine",
						 				 "Configuration profile": "$workflow.profile",
								 		 "Start timestamp": "$workflow.start",
		 								 "Stop timestamp": "$workflow.complete",
										 "Run name": "$workflow.runName",
										 "Revision": "$workflow.revision",
										 "Commit ID": "$workflow.commitId",
										 "Duration": "$workflow.duration"]]
	def params_data = ["Parameters": [:]]
	for (p in params) {
		if (!p.key.contains('-')) {
			params_data."Parameters".put("$p.key", "$p.value")
		}
		if (p.key.equals('testing_run')) {
			testing = true
			data."Workflow Information".put("Exit status", "$workflow.exitStatus")
			data."Workflow Information".put("Error message", "$workflow.errorMessage")

		}
		
	}

	def manifest_data = ["Manifest": [:]]
	for (p in workflow.manifest.getProperties()) {
		p = p.toString()
		def split_str = p.split("=")
		if (split_str[0].equals("name") || split_str[0].equals("version")) {
			manifest_data."Manifest".put(split_str[0], split_str[1])
		}
	}

	def json_str = JsonOutput.toJson(data+params_data+manifest_data)
	def json_beauty = JsonOutput.prettyPrint(json_str)
	workflow_info = file("$params.outDir/workflow_info.json")
	workflow_info.write(json_beauty)
}
