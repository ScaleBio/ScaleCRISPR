def shouldRunLibDetection(samples, quantum) {
	def header = samples.withReader { reader ->
		reader.readLine()
	}
	if (!header.contains("libName") && quantum && params.fastqDir) {
		return true
	}
	return false
}

// Create a 'file()' from string 'path'
// 'path' could be null, a URL, an absolute file-path or relative to 'baseDir'
def expandPath(path, baseDir) {
	if (path == null) { return null}
	def uri = new URI(path)
	if (uri.isAbsolute()) {
		return file(path)
	} else {
		return baseDir.resolve(path)
	}
}

def resolveCrisprToRnaMapping(libJson) {
	if (params.crisprToRnaMapping) {
		crisprToRnaMapping = file(params.crisprToRnaMapping, checkIfExists:true)
	} else {
		libStructure = Utils.loadJson(libJson)
		if (libStructure.get("crispr_to_rna_mapping")) {
			crisprToRnaMapping = file(libJson.getParent()) / libStructure["crispr_to_rna_mapping"]
		} else {
			// path input cannot be null, no mapping file used in 3lvl pipeline
			crisprToRnaMapping = []
		}
	}
	return crisprToRnaMapping
}