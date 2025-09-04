import nextflow.Nextflow
import java.net.URI;


// Adapted from nf-core/rnaseq
class ParamLogger {
    public static void initialise(workflow, params, log) {
        if (params.help) {
            log.info("\nScaleBio Seq Suite: RNA Workflow\nFor usage information, please see README.md\n")
            System.exit(0)
        }
        computeOutdir(workflow, params)
        validateParams(params, log)
    }

    public static void computeOutdir(workflow, params) {
        def uriOutDir
        if (params.computeOutDir) {
            // Add trailing slash before creating a URI so that is interpreted as a directory and not a file
            // This is necessary for URI.resolve to work correctly
            // If a trailing slash already exists, adding an extra one does not create an issue
            uriOutDir = URI.create(params.outDir + "/")
            if (params.runFolder) {
                def uriRunFolderName = params.runFolder.tokenize("/").last()
                params.outputDir = uriOutDir.resolve(uriRunFolderName+"/").resolve(workflow.runName+"/").toString()
            }
            if (params.fastqDir) {
                def uriFastqDirName = params.fastqDir.tokenize("/").last()
                params.outputDir = uriOutDir.resolve(uriFastqDirName+"/").resolve(workflow.runName+"/").toString()
            }
            if (params.reporting) {
                if (params.resultDir) {
                    def uriResultDirName = params.resultDir.tokenize("/").last()
                    params.outputDir = uriOutDir.resolve(uriResultDirName+"/").resolve(workflow.runName+"/").toString()
                } else {
                    params.outputDir = uriOutDir.resolve(workflow.runName+"/").toString()
                }
            }
            if (params.ultimaCramDir) {
                def uriUltimaCramDirName = params.ultimaCramDir.tokenize("/").last()
                params.outputDir = uriOutDir.resolve(uriUltimaCramDirName+"/").resolve(workflow.runName+"/").toString()
            }
        } else {
            params.outputDir = params.outDir
        }
    }

    public static String throwError(errMessage) {
        Map colors = logColours()
        Nextflow.error("${colors.red}ERROR${colors.reset}: $errMessage")
    }
    
    // Takes a camelCase string and converts it to kebab-case
    // Necessary for validating parameters since nextflow converts camelCase to kebab-case
    public static String camelToKebab(String camelCase) {
        return camelCase.replaceAll(/([a-z])([A-Z])/, '$1-$2').toLowerCase()
    }

    public static checkArg(arg, str) {
        if (arg == "") {
            throwError("Parameter $str has invalid value of ''")
        }
        if (arg == true) {
            throwError("Parameter $str has invalid value of true")
        }
        if (arg == false) {
            throwError("Parameter $str has invalid value of false")
        }
    }

    // Check required workflow inputs
    public static validateParams(params, log) {
        def parametersToCheck = ['samples', 'runFolder', 'fastqDir', 'libStructure']
        parametersToCheck.each { param ->
            checkArg(params.getProperty(param), param)
        }
        if (params.samples == null || params.samples == true) {
            throwError("Must specify --samples <samples.csv>")
        }
        if (params.fastqDir && params.runFolder) {
            throwError("Cannot specify both --runFolder and --fastqDir")
        }
        if (params.libStructure == null || params.libStructure == true) {
            throwError("Must specify --libStructure <libQuantumV1.0.json>")
        }
        if (params.reporting) {
            if (params.runFolder || params.fastqDir) {
                throwError("Cannot specify --runFolder or --fastqDir when running reporting-only (--reporting)")
            }
        }
    }


    public static Map logColours(monochrome_logs=false) {
        Map colorcodes = [:]
        // Reset / Meta
        colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
        colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
        colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
        colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
        colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
        colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
        colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"

        // Regular Colors
        colorcodes['black']      = monochrome_logs ? '' : "\033[0;30m"
        colorcodes['red']        = monochrome_logs ? '' : "\033[0;31m"
        colorcodes['green']      = monochrome_logs ? '' : "\033[0;32m"
        colorcodes['yellow']     = monochrome_logs ? '' : "\033[0;33m"
        colorcodes['blue']       = monochrome_logs ? '' : "\033[0;34m"
        colorcodes['purple']     = monochrome_logs ? '' : "\033[0;35m"
        colorcodes['cyan']       = monochrome_logs ? '' : "\033[0;36m"
        colorcodes['white']      = monochrome_logs ? '' : "\033[0;37m"

        // Bold
        colorcodes['bblack']     = monochrome_logs ? '' : "\033[1;30m"
        colorcodes['bred']       = monochrome_logs ? '' : "\033[1;31m"
        colorcodes['bgreen']     = monochrome_logs ? '' : "\033[1;32m"
        colorcodes['byellow']    = monochrome_logs ? '' : "\033[1;33m"
        colorcodes['bblue']      = monochrome_logs ? '' : "\033[1;34m"
        colorcodes['bpurple']    = monochrome_logs ? '' : "\033[1;35m"
        colorcodes['bcyan']      = monochrome_logs ? '' : "\033[1;36m"
        colorcodes['bwhite']     = monochrome_logs ? '' : "\033[1;37m"

        // Underline
        colorcodes['ublack']     = monochrome_logs ? '' : "\033[4;30m"
        colorcodes['ured']       = monochrome_logs ? '' : "\033[4;31m"
        colorcodes['ugreen']     = monochrome_logs ? '' : "\033[4;32m"
        colorcodes['uyellow']    = monochrome_logs ? '' : "\033[4;33m"
        colorcodes['ublue']      = monochrome_logs ? '' : "\033[4;34m"
        colorcodes['upurple']    = monochrome_logs ? '' : "\033[4;35m"
        colorcodes['ucyan']      = monochrome_logs ? '' : "\033[4;36m"
        colorcodes['uwhite']     = monochrome_logs ? '' : "\033[4;37m"

        // High Intensity
        colorcodes['iblack']     = monochrome_logs ? '' : "\033[0;90m"
        colorcodes['ired']       = monochrome_logs ? '' : "\033[0;91m"
        colorcodes['igreen']     = monochrome_logs ? '' : "\033[0;92m"
        colorcodes['iyellow']    = monochrome_logs ? '' : "\033[0;93m"
        colorcodes['iblue']      = monochrome_logs ? '' : "\033[0;94m"
        colorcodes['ipurple']    = monochrome_logs ? '' : "\033[0;95m"
        colorcodes['icyan']      = monochrome_logs ? '' : "\033[0;96m"
        colorcodes['iwhite']     = monochrome_logs ? '' : "\033[0;97m"

        // Bold High Intensity
        colorcodes['biblack']    = monochrome_logs ? '' : "\033[1;90m"
        colorcodes['bired']      = monochrome_logs ? '' : "\033[1;91m"
        colorcodes['bigreen']    = monochrome_logs ? '' : "\033[1;92m"
        colorcodes['biyellow']   = monochrome_logs ? '' : "\033[1;93m"
        colorcodes['biblue']     = monochrome_logs ? '' : "\033[1;94m"
        colorcodes['bipurple']   = monochrome_logs ? '' : "\033[1;95m"
        colorcodes['bicyan']     = monochrome_logs ? '' : "\033[1;96m"
        colorcodes['biwhite']    = monochrome_logs ? '' : "\033[1;97m"

        return colorcodes
    }
}