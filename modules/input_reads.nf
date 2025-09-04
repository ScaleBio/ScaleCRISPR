/*
* Perform fastq generation, barcode processing, sample demux and qc of input data
*
*/
include { LibraryDetection } from './library_detection.nf'

def getReadFromFqname(fqName) {
    def m = ( fqName =~ /_([RI][12])[_.]/ )
    if (!m) { return null }
    return m[0][1]
}

// Return the libraries that were filtered out
def getFilteredLibraries(allLibNames, allInputLibNames, reason) {
	allInputLibNames
		.map { it ->
			def rows = []
			for (libName in it) {
				// Indicates library that was filtered out due to @reason
				if (!(libName in allLibNames.val)) {
					rows.add(",${libName},0,${reason}")
				}
			}
			rows
		}
}

// Extract sampleName (and optionally subsample, split by RT well) for one bcParser output
// @bcParserOut -> library name, fastq file
def constructSampleName (bcParserOut, splitFastq) {
	def libName = bcParserOut[0]
	def file = bcParserOut[1]
	def tok = file.getName().toString().tokenize('_') // bcParser outputs are "sampleName_well_...fastq.gz"
	def sampleName = tok[0]
	def subsampleName = sampleName
	if (splitFastq) {
		subsampleName = "${sampleName}_${tok[1]}"
	}
	return tuple(libName, sampleName, subsampleName, file)
}

def constructMatchingKey(fname) {
	def identifying_char_removed = fname.replaceAll(/_[RI][12][_\.]/, "")
	return identifying_char_removed
}

// Create a bcl-convert samplesheet for libraries in samples.json
process MakeBclConvertSamplesheet {
	label 'small'

	publishDir file(params.outDir) / "fastq", mode: 'copy'

	input: 
	path(samplesCsv)
	path(libStructDir) // Directory containing the library structure definition (barcode sequence lists are loaded from here)
	val(libStructName) // Filename of the library structure definition json
	path(runinfo) // RunInfo.xml from Sequencer RunFolder (Read-lengths etc.)

	output: 
	path("samplesheet.csv")

	script:
	libStruct = "$libStructDir/$libStructName"
	opts = ""
	if (params.splitFastq) {
		opts += "--splitFastq "
	}
	"""
	bcl_convert_sheet.py $samplesCsv $libStruct $runinfo $opts > samplesheet.csv
	"""
}

/*Run bcl-convert, used when starting from a sequencing run-folder
  Requires a separate bcl-convert samplesheet*/
process BclConvert {
	publishDir params.outDir, pattern: 'fastq/Reports/*', mode: 'copy'
	publishDir params.outDir, pattern: 'fastq/*.fastq.gz', enabled: params.fastqOut

	input: 
	path(run)
	path(samplesheet)

	output: 
	path("fastq/*fastq.gz"), emit: fastq
	path("fastq/Reports/*"), emit: stats

	script:
	wthreads = (task.cpus/2 - 0.5).round()
	dthreads = (task.cpus/3 - 0.5).round()
	"""
	bcl-convert --sample-sheet $samplesheet --bcl-input-directory $run \
    --bcl-num-conversion-threads $wthreads --bcl-num-compression-threads $wthreads --bcl-num-decompression-threads $dthreads \
    $params.bclConvertParams --output-directory fastq
	"""
}

// Run fastQC on (a pair of) fastq files from inputs or bcl-convert
// Note that these are the input fastq files, not outputs of bcParser
process FastQC {
	tag "${fq[0].getSimpleName()}"

	publishDir file(params.outDir) / "fastq", mode: 'copy'

	input:
	path(fq)
	
	output:
	path("fastqc/*.html"), emit: html
	path("fastqc/*.zip"), emit: zip

	"""
	mkdir fastqc
	fastqc -t 2 -o fastqc $fq
	"""
}

// Use multiQC to report a single QC report for all fastQC and bcl-convert reports
process MultiQC {
	label 'optional'
	memory { (reports.size() > 100 ? 12.GB : 6.GB) * task.attempt}

	publishDir file(params.outDir) / "reports", mode: 'copy'
	
	input:
	path(reports)
	
	output:
	path("multiqc_report.html")

	"""
	multiqc .
	"""
}

// Run bcParser to extract and correct cell-barcodes
// Optionally demuxes fastq files and produces a bam file based on some barcodes
process BarcodeParser {
	tag "$libName ${task.index}"

	publishDir { file(params.outDir) / dir_name }, pattern: "$demuxDir/*bam", enabled: params.bcParserBamOut
	publishDir { file(params.outDir) / dir_name }, mode: 'copy', pattern: "$demuxDir/*{txt,tsv,json}"
	/* Publish barcode demuxed bam files if flag is turned on; Always publish metrics files */

	input:
	// samples.csv
	path(sheet)
	
	// Directory containing the library structure definition (barcode sequence lists are loaded from here)
	path(libStructDir)
	// Filename of the library structure definition json
	val(libStructName)
	// Input fastq file
	tuple val(libName), path(fqFiles)
	// Fastq files are reverse complemented or not
	val(isRc)

	output:
	//Output bam might not exist if the input has no passing reads
	tuple val(libName), path("$demuxDir/*_S[1-9]*bam"), emit: bam, optional: true
	path("$demuxDir/*_S0*bam"), emit: unknown, optional: true
	path("$demuxDir/*.tsv")
	tuple(val(libName), path("$demuxDir/metrics.json"), emit: metrics)

	script:
	opts = ""
	demuxDir = "${libName}.demux"
	libStruct = "$libStructDir/$libStructName"
	dir_name = "barcodes"
	if (params.splitFastq) {
		dir_name = "$dir_name/split_bcparser_jobs/bcparser.${libName}.${task.index}"
	}
	if (isRc) {
		opts += "--rev-comp-idx2 "
	}
	"""
	bc_parser --lib-struct $libStruct --demux $sheet --lib-name $libName -v --reads ${fqFiles.join(" ")} --out $demuxDir --write-bam $opts
	"""
}

process MergeDemux {
	tag "$libName"
	label 'small'
	
	publishDir { file(params.outDir) / "barcodes" }, mode:'copy'

	input:
   	tuple(val(libName), path("demux_metrics*.json"))
	path(libStructDir)
	val(libStructName)

	output:
    tuple(val(libName), path("*metrics.json"))

	script:
	libStruct = "$libStructDir/$libStructName"
	"""
    merge_bc_parser_output.py --bc_jsons demux_metrics* --lib_json $libStruct --libName $libName
    """
}

process SeqTagFinder {
	tag "$subsample"

	input:
	tuple val(sample), val(subsample), path(untaggedBam), path(guideWhitelist)

	output:
	tuple val(sample), val(subsample), path("$dir_name/*_S[1-9]*bam"), emit: bam, optional: true

	script:
	dir_name = "taggedBams"
	"""
	seq_tag_finder --bams ${untaggedBam.join(" ")} --whitelist $guideWhitelist --tag_in_output_bam tg
	"""
}

// Fastq generation, trimming, QC, barcode extraction and sample demux
workflow INPUT_READS {
take:
	samples
	samplesCsv
	libraryInfo
	runFolder
	fastqDir
	runLibDetection
main:
	runDir = runFolder
	fqDir = fastqDir
	fqs = null
	samples
		.collect { it.libName }
		.unique()
		.set { allInputLibNames }
	if (runDir != null) {
		if (params.fastqSamplesheet == null) {
			MakeBclConvertSamplesheet(
				samplesCsv,
				libraryInfo.rnaLibraryStructureFile.getParent(),
				libraryInfo.rnaLibraryStructureFile.getName(),
				file("$runDir/RunInfo.xml", checkIfExists:true)
			)
			fqSheet = MakeBclConvertSamplesheet.out
		} else {
			fqSheet = file(params.fastqSamplesheet)
		}
		BclConvert(
			file(runDir),
			fqSheet
		)
		fqs = BclConvert.out.fastq.flatten()
	} else if (fqDir != null) {
		fqs = Channel.fromPath("$fqDir/**fastq.gz", checkIfExists: true)
	} else {
		ParamLogger.throwError("Must specify either '--runFolder' or '--fastqDir' when running alignment")
	}
	// fqs -> (fastq file)
	fqs.dump(tag:'fqs')

	def isRc = false
	// Need to set empty channel for the 3lvl case
	filteredLibs = channel.value([])
	if (runLibDetection) {
		fqs
			.filter{ !it.getName().contains("Undetermined") }
			.filter{ getReadFromFqname(it.getName()) == 'I2' } // file containing the library barcode
			.ifEmpty { ParamLogger.throwError("No index2 fastq files found") }
			.buffer(size: params.filesPerLibDetection, remainder: true)
			.dump(tag:'index2Fq')
			.set { idx2Fq }
		// Run library detection on all index2 files
		LibraryDetection(
			idx2Fq,
			libraryInfo.rnaLibraryStructureFile.getParent(),
			libraryInfo.rnaLibraryStructureFile.getName(),
		)
		LibraryDetection.out.pool_mapping
			.collectFile(
				name: 'fastq_file_to_library_assignment.csv',
				keepHeader: true,
				storeDir: file(params.outDir) / "fastq"
        	)
			.splitCsv(header:true, strip:true)
			.filter { row -> row.pool != "Unknown" && row.pool != "Ambiguous" }
			.collect()
			.map { rows ->
				def is_rc_values = rows.collect { it.is_rc }.unique()
				if (is_rc_values.size() > 1) {
					ParamLogger.throwError("Input fastq files should all be either reverse complemented or not")
				}
				rows // Return all rows if check passes
			}
			.flatMap { it } // Flatten the list of rows back to individual rows
			.dump(tag:'poolMapping')
			.set { poolMapping }
		poolMapping
			// Can look at the is_rc value of any row since all rows will have the same value
			.first()
			.filter { !it.is_rc }
			.ifEmpty { isRc = true }
		poolMapping
			.map { it ->
				// Remove I2 from fastq file name to enable joining with other fastq files of same library
				tuple(constructMatchingKey(it.file), it.pool)
			}
			// Cross because for one I2 file, there is one I1, one R1 and one R2 file 
			.cross(
				// fqs corresponds to all files, whereas poolMapping just corresponds to I2
				fqs
					.map { file ->
						tuple(constructMatchingKey(file.getName().toString()), file)
					}
			)
			.map { it -> //[[matching_key, pool], [matching_key, file]]
				tuple(it[0][1], it[0][0], it[1][1])
			}
			.set { fqFiles }
		poolMapping
			.map { it.pool }
			.unique()
			.collect()
			.ifEmpty { ParamLogger.throwError("No valid PCR pool found in index2 fastq files") }
			.dump(tag:'validPools')
			.set { validPools }
		samples
			// Filter out samples that have a library which was not found during library detection
			.filter { sample ->
				validPools.val.contains(sample.libIndex2)
			}
			.dump(tag:'samplesAfterLibraryDetection')
			.set { samples }
		samples
			.collect { it.libName }
			.unique()
			.set { allLibNamesAfterLibraryDetection }
		filteredLibs = getFilteredLibraries(allLibNamesAfterLibraryDetection, allInputLibNames, "Library not found in input fastq files")
	} else {
		if (libraryInfo.quantum) {
			fqs
				.filter{ getReadFromFqname(it.getName()) == 'I2' }
				.ifEmpty { ParamLogger.throwError("No index2 fastq files found") }
				.filter { it.size() < params.index2MinFileSize } // Get files that are less than 1MB to filter out later
				.collect()
				.map { it ->
					def matchingKeyForReadsLessThanOneMB = []
					for (file in it) {
						matchingKeyForReadsLessThanOneMB << constructMatchingKey(file.getName().toString())
					}
					matchingKeyForReadsLessThanOneMB
				}
				.dump(tag:'fqsToFilter')
				.set { fqsToFilter }
			fqsToFilter.collect().view { items ->
				log.info "Filtering out fastq files ${items} that are less than $params.index2MinFileSize bytes"
			}
			fqs
				.collect()
				.map { allFastqs ->
					def passingFastqs = [] 
					for (fastq in allFastqs) {
						def matchingKey = constructMatchingKey(fastq.getName().toString())
						if (!(matchingKey in fqsToFilter.val)) {
							passingFastqs << fastq
						}
					}
					passingFastqs
				}
				.flatMap()
				.dump(tag:'passingFastqsPostFileSizeFilter')
				.set { fqs }
			fqs
				.map { file ->
					file.getName().toString().tokenize('_')[0]
				}
				.collect()
				.set { allLibNames }
			samples
				.filter { it -> allLibNames.val.contains(it.libName) }
				.ifEmpty { ParamLogger.throwError("No matching fastq files found for libraries") }
				.dump(tag:'samplesAfterFilteringSmallFiles')
				.set { samples }
			samples
				.collect { it.libName }
				.unique()
				.set { allLibNamesAfterFilteringSmallFiles }
			filteredLibs = getFilteredLibraries(
				allLibNamesAfterFilteringSmallFiles,
				allInputLibNames,
				"Library filtered because all index2 fastq files are less than $params.index2MinFileSize bytes"
			)
		}
		// Organize fastq files by sample
		// fqFiles_with_fastq_matching_key -> (library_name, fastq_matching_key, file)
		fqFiles = fqs.map { file ->
			def fname = file.getName().toString()
			def libName = fname.tokenize('_')[0]
			// Construct matching key to join on later on in the workflow
			if (params.splitFastq) {
				return tuple(libName, constructMatchingKey(fname), file)
			}
			// Matching key is just the library name if params.splitFastq is false
			else {
				return tuple(libName, libName, file)
			}
		}
	}
	fqFiles.dump(tag:'fqFiles')
	pairedFqFiles = fqFiles.groupTuple(by:[0,1])
	// This keeps only those sets of fastq files that correspond to a sample from our samples.csv
	// it.libName takes only the library name from samples csv
	// unique because samples csv has multiple entries for multiple samples, but only one library name
	// cross ensures that the right library names are paired with the right fastq files
	// The last map is to remove the extra library name element(it[0]) in the tuple that we get after the cross
	fqSamples = samples.map({it.libName}).unique().cross(pairedFqFiles).map{it[1]}
	// fqSamples -> (library name, [I1, R1, R2])
	fqSamples.dump(tag:'fqSamples')

	// Check that for each library (unique libName in samples.csv), we have a complete set of fastq files
	// checkFastq -> (libName, sampleName, Fastqs)
	// This slightly strange channel construction (including it.id) works around the lack of 'left-join'
        // in nextflow (this way we can distinguish 'left-only' and 'right-only' join results)
	checkFastq = samples.map{ [it.libName, it.id] }.unique{ it[0] }.join(fqSamples, remainder: true)
	checkFastq.dump(tag:'checkFastq')
	checkFastq.map {
		// If samples.csv has a library name which does not have corresponding fastq files
		// then 3rd element of checkFastq will be null
		if (it[3] == null) {
			ParamLogger.throwError("Library ${it[0]} does not have matching fastq files. None of the provided fastq filenames start with ${it[0]}")
		}
		// Check that each library has index1, read1 and read2 fastq files
		if (it[3].any {fq -> fq.getName().contains("_R1_") || fq.getName().contains("_R1.")} == false) {
			ParamLogger.throwError("Library ${it[0]} does not have read1 fastq file")
		}
		if (it[3].any {fq -> fq.getName().contains("_R2_") || fq.getName().contains("_R2.")} == false) {
			ParamLogger.throwError("Library ${it[0]} does not have read2 fastq file")
		}
		if (it[3].any {fq -> fq.getName().contains("_I1_") || fq.getName().contains("_I1.")} == false) {
			ParamLogger.throwError("Library ${it[0]} does not have index fastq file")
		}
	}

	// Get rna lib names
	samples
		.collect { it.libName }
		.dump(tag:'rnaLibNames')
		.set { rnaLibNames }

	// Get all remaining non hash fastq files
	rnaLibNames
		.flatMap()
		.unique()
		.cross(
			fqSamples
				.map {
					def read = []
					// Iterate over all fastq files associated with a libName
					for (fq in it[2]) {
						read << fq
					}
					[it[0], it[1], read]
				}
		)
		.map { it[1] } // TODO drop matching key here
		// fqSamplesRnaReads -> (library name, matching key, [R1, R2, I1, I2])
		.dump(tag:'fqSamplesRnaReads')
		.set { fqSamplesRnaReads }
	
	// Counter that we use for grouping bcParser jobs with the same libName together
	// Intention is to create params.bcParserJobs number of groups of fastq files belonging to the same libName
	// count%params.bcParserJobs will give us a number between 0 and params.bcParserJobs-1 which we can then group on
	def count = 0
	// Join the trimmed fastq files with hash fastq files and rna fastq files
	fqSamplesRnaReads
		// Sort channel so results are cached
		.toSortedList { a, b ->
			// Sort by libName and then fastq file names
			def result = a[0] <=> b[0]
			if (result == 0) {
				result = a[2].join(",") <=> b[2].join(",")
			} 
			return result
		}
		.flatMap()
		.map {
			tuple(it[0], it[2], ++count % params.bcParserJobsPerLibName)
		}
		// Group on libName and counter
		.groupTuple(by:[0,2], sort: { a, b -> a.join(",") <=> b.join(",") })

		// Concatenate list of lists into single list
		.map { lib, fastqFiles, _count ->
			[lib, fastqFiles.flatten()]
		}
		.dump(tag:'bcParserInput')
		.set { bcParserInput }

	// Process cell-barcodes and (optionally) split fastqs into samples based on sample barcode
	BarcodeParser(
		samplesCsv,
		libraryInfo.rnaLibraryStructureFile.getParent(),
		libraryInfo.rnaLibraryStructureFile.getName(),
		bcParserInput,
		isRc
	)
	BarcodeParser.out.bam.dump(tag:'bcParserBam')
	// Check that there exists bam files post bcParser, for every sample that corresponds to a libName in samples.csv
	samples
		.map{ [it.libName, it.id] }
		// get all samples corresponding to a libName
		.groupTuple()
		.join(BarcodeParser.out.bam.groupTuple(), remainder: true)
		// bcParserOutCheck -> (libName, [sample names corresponding to that libName], [[output bam files]])
		.dump(tag:'bcParserOutCheck')
		.map {
			if (it[2] == null) {
				ParamLogger.throwError("Library ${it[0]} does not have any passing bam files after bcParser")
			}
			// iterate over every sample name corresponding to a libName
			for (sample in it[1]) {
				// extract the sample name from the id
				def sampleName = sample.tokenize('.')[0]
				def exists = false
				// iterate over every bam file for that libName and check that the sample name is present in the bam file name
				for (file in it[2].flatten()) {
					if (file.getName().toString().contains(sampleName)) {
						exists = true
					}
				}
				if (!exists) {
					log.warn("Sample ${sampleName} does not have passing bam files from bcParser.")
				}
			}
		}
	// Calculate total sample reads for rna and scaleplex libraries
	rnaLibNames
		.flatMap()
		.unique()
		.cross(BarcodeParser.out.metrics)
		.collect()
		.map { it -> // it -> [libName1, [libName1, metrics], libName2, [libName2, metrics], ...]
			def readsBySampleID = [:]
			for (int i = 0; i < it.size(); i += 2) {
				def libName = it[i]
				// Don't need to check for i+1 being out of bounds since we know that it.size() is a multiple of 3
				def metricsJson = Utils.loadJson(it[i+1][1])
				for (sample in metricsJson['samples']) {
					// Groovy returns the entire map, along with the key when using the nomenclature above (for key in map)
					// value returns the contents of the map
					def sampleMap = sample.value
					if (sampleMap['name'] == "Unknown") {
						continue
					}
					def totalReads = sampleMap['passingreads'].toLong() + sampleMap['tooshort'].toLong()
					def id = sampleMap['name'] + "." + libName
					if (id in readsBySampleID) {
						readsBySampleID[id]['totalReads'] += totalReads
						readsBySampleID[id]['passingreads'] += sampleMap['passingreads'].toLong()
					} else {
						readsBySampleID[id] = ['totalReads': totalReads, 'passingreads': sampleMap['passingreads'].toLong()]
					}
				}
			}
			def outputFile = file("${params.outDir}/reports/reads_per_sample.csv")
			// Parent directory does not exist at this point
			outputFile.getParent().mkdirs()
			// Sort so that csv has rows in deterministic order
			def sortedReadsBySampleID = readsBySampleID.entrySet().sort { a, b -> a.key <=> b.key }
			def contentToWrite = ["Sample,Library,PassingReads,Comment"]
			for (entry in sortedReadsBySampleID) {
				def sampleName = entry.key.tokenize('.')[0]
				def libName = entry.key.tokenize('.')[1]
				def msg = ""
				if ((libName in rnaLibNames.val) && (entry.value['passingreads'] < params.minPassingSampleReads)) {
					msg = "Filtered because sample has less than ${params.minPassingSampleReads} passing reads"
				}
				contentToWrite.add("${sampleName},${libName},${entry.value['passingreads']},${msg}")
			}
			outputFile.text = (contentToWrite + filteredLibs.val).join('\n')
			// Drop samples where the passingreads is less than minPassingSampleReads
			readsBySampleID = readsBySampleID.findAll { key, value -> value['passingreads'] >= params.minPassingSampleReads }
			def totalSampleReadsBySampleID = [:]
			// Need only totalReads for computation downstream, drop passingreads
			for (sample in readsBySampleID) {
				totalSampleReadsBySampleID[sample.key] = sample.value['totalReads']
			}
			// Convert map to list of lists to enable converting it to a channel that can be joined on later
			totalSampleReadsBySampleID.collect { key, value -> [key, value] }
		}
		.flatMap()
		.dump(tag:'totalSampleReadsBySampleID')
		.set { totalSampleReadsBySampleID }
	// Get all sample IDs that have passing reads
	totalSampleReadsBySampleID
		.map { it[0] }
		.collect()
		.ifEmpty { ParamLogger.throwError("No RNA sample in experiment has ${params.minPassingSampleReads} reads post barcode demux") }
		.set { sampleIdFromTotalSampleReadsBySampleID }
	// Construct samples object with samples that have passing reads
	samples
		.filter { it -> sampleIdFromTotalSampleReadsBySampleID.val.contains(it.id) }
		.dump(tag:'samplesWithPassingReads')
		.set { samplesWithPassingReads }
	// Get libNames from samples that have passing reads
	samplesWithPassingReads
		.collect { it.libName }
		.dump(tag:'rnaLibNamesWithNonZeroReads')
		.set { rnaLibNamesWithNonZeroReads }
	
	BarcodeParser.out.bam
		.transpose()
		.map{ constructSampleName(it, params.splitFastq) }
		.dump(tag:'bcBam')
		.set { bcBam }

	bcBam
		.combine(rnaLibNamesWithNonZeroReads.flatMap().unique(), by:0) // only rna reads
		.combine(samplesWithPassingReads.map{ [it.libName, it.sample, it.id, file(params.guidesDir) / it.guide] }, by:[0,1])
		// id, subsample, bam, guide
		.map{ [it[4], it[2], it[3], it[5]] }
		.dump(tag:'bcBamRna')
		.set { bcBamRna }

	if (params.fastqc) {
		fqSamples
			.flatMap{it.get(2)}
			// Exclude I1 and I2 from fastqc
			.filter{getReadFromFqname(it.getName())[0] == 'R'}
			// Batch R1 and R2 into one fastqc job
			.map { fastq_file ->
				def fname = fastq_file.getName().toString()
				// Get key for grouping R1 and R2 together
				def key = fname.replaceAll(/_R[12]/, '')
				tuple(key, fastq_file)
			}
			.groupTuple(size: 2, remainder: true)
			.map { _key, fastq_files ->
				// Check that we have exactly 2 fastq files for each key
				if (fastq_files.flatten().size() != 2) {
					ParamLogger.throwError("Fastq files for ${_key} are not paired")
				}
				fastq_files.flatten()
			}
			.collect()
			.map { it ->
				// If number of R1, R2 pairs is greater than totalFastqcJobs, split them into smaller groups
				if (it.size() / 2 > params.totalFastqcJobs) {
					// Math.ceil to ensure number of jobs stays within params.totalFastqcJobs
					def group_size = (int) Math.ceil(it.size() / params.totalFastqcJobs)
					if (group_size % 2 == 1) {
						group_size += 1 // Ensure group size is even for R1 and R2 pairing
					}
					return it.collate(group_size)
				}
				// For the case where we don't need to do any more batching
				it.collate(2)
			}
			.flatMap()
			.dump(tag:'mainReads')
			.set{ mainReads }
		FastQC(mainReads)
		reports = FastQC.out.zip
		if (runDir != null) {
			reports = reports.mix(BclConvert.out.stats)
		}
		MultiQC(reports.collect())
	}

	if (params.splitFastq) {
		// Need to merge the bcParser outputs for a single library together
		organized_demux = BarcodeParser.out.metrics.groupTuple()
		MergeDemux(
			organized_demux,
			libraryInfo.rnaLibraryStructureFile.getParent(),
			libraryInfo.rnaLibraryStructureFile.getName(),
		)
		metrics = MergeDemux.out
    }
	else {
		metrics =  BarcodeParser.out.metrics
	}
	SeqTagFinder(bcBamRna)
	SeqTagFinder.out.bam.dump(tag:'seqTagFinderBam')
emit:
	ubam = SeqTagFinder.out.bam
	metrics = metrics
	totalSampleReadsBySampleID = totalSampleReadsBySampleID
	samples = samplesWithPassingReads
}
