//figure generation
process CropAnalysis {
tag "$sample"
label 'crispr'
publishDir { file(params.outDir) / "reports" }, mode: 'copy'
input:
	val(libJsonName)
	path(qc_template)
	tuple val(sample), path(guideTab), path(cellMetrics), path("metrics.csv"), val(barcodesRange)
output:
	path("plots/*.png")
	path("csv/*.txt")
	path("${sample}.report.html")
script:
"""
	cropAnalysis.py --sample ${sample} --umi_matrix ${guideTab} --guideMetricsCsv ${cellMetrics} \
		--thresh ${params.thresh} --outDir ./ --template ${qc_template} --libJsonName $libJsonName \
		--barcodesRange '$barcodesRange' --workflowVersion $workflow.manifest.version
"""
}

process SampleStats {
tag "$sample"
publishDir { file(params.outDir) / "samples" / sample }, pattern: "metrics.csv", mode: 'copy', saveAs: { "${sample}.sample.metrics.csv" }

input:
tuple val(sample), path("cellMetrics_???.csv"), path("cellMetricsFiltered_???.csv"), path(guides)

output:
tuple val(sample), path("metrics.csv"), emit: sample_metrics

script:
"""
	sample_stats.py --sample $sample --cellMetrics cellMetrics_*.csv --cellMetricsFiltered cellMetricsFiltered_*.csv \
		--numGuides \$(wc -l < $guides) --thresh ${params.thresh} --threads $task.cpus --memory '$task.memory'
"""
}