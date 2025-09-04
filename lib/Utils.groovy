class Utils {
    // Load .json file into object
    // (reference genomes, library design, etc.)
    public static Map loadJson(json) {
        def jsonSlurper  = new groovy.json.JsonSlurper()
        return jsonSlurper.parse( json )
    }

    // Derive value for the 'barcodes' field in the merged sample report from 
    // a list of the 'barcodes' column value for each sub-sample
    // Usually all samples have the same barcode entry (multi-plate use-case)
    public static String combineSampleBarcodes(barcodes) {
        ArrayList bcs = barcodes.unique()
        if (bcs.size() == 1) {
            return bcs[0]
        } else {
            return "combined"
        }
    }
}
