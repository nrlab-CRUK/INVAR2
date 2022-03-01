@Grab('com.xlson.groovycsv:groovycsv:1.3')

/**
 * Validation pipeline configuration.
 * Make sure mandatory parameters are defined and all others are within acceptable ranges.
 */

import java.util.regex.Pattern
import static java.nio.file.Files.exists
import static org.apache.commons.lang3.StringUtils.isBlank
import static org.apache.commons.lang3.StringUtils.trimToEmpty

/**
 * Check if value is a boolean.
 */
def isBoolean(value)
{
    return value instanceof Boolean
}

/**
 * Check if value is acceptable as an integer.
 */
def isInteger(value)
{
    return value instanceof Short || value instanceof Integer || value instanceof Long
}

/**
 * Check if value is acceptable as a decimal.
 * Any number is acceptable for this, including integers.
 */
def isNumber(value)
{
    return value instanceof Number
}

/**
 * Find a file on a Unix/Windows style path list. The path is
 * searched in order, and the first file found with an exact name
 * match is returned. Returns null if there is no file with the name
 * on any path.
 */
def findFileOnPath(path, filename)
{
    def pathParts = path.split(Pattern.quote(File.pathSeparator))
    def dirFiles = []

    for (dirPath in pathParts)
    {
        def dir = new File(dirPath)
        if (dir.directory)
        {
            dirFiles << dir
        }
    }

    for (dir in dirFiles)
    {
        def file = new File(dir, filename)
        if (file.exists())
        {
            return file
        }
    }

    return null
}


def validateParameters(params)
{
    def errors = false

    params.with
    {
        if (!containsKey('STUDY'))
        {
            log.error "STUDY must be set in your project's configuration file."
            errors = true
        }
        if (!containsKey('ERROR_SUPPRESSION_NAME'))
        {
            log.error "ERROR_SUPPRESSION_NAME must be set in your project's configuration file."
            errors = true
        }

        if (!containsKey('FAMILY_SIZE'))
        {
            log.error "FAMILY_SIZE must be set in your project's configuration file."
            errors = true
        }
        else
        {
            if (isInteger(FAMILY_SIZE))
            {
                if (FAMILY_SIZE < 1)
                {
                    log.error "FAMILY_SIZE must be a positive integer > 0."
                    errors = true
                }
            }
            else
            {
                log.error "FAMILY_SIZE must be an integer."
                errors = true
            }
        }

        /**
         * Input files.
         */

        if (!containsKey('LAYOUT_TABLE'))
        {
            log.error "LAYOUT_TABLE must be set in your project's configuration file."
            errors = true
        }
        else
        {
            def layoutFile = file(LAYOUT_TABLE)
            if (!exists(layoutFile))
            {
                log.error "LAYOUT_TABLE file ${layoutFile.toAbsolutePath()} does not exist."
                errors = true
            }
        }
        if (!containsKey('TUMOUR_MUTATIONS_CSV'))
        {
            log.error "TUMOUR_MUTATIONS_CSV must be set in your project's configuration file."
            errors = true
        }
        else
        {
            def mutationsFile = file(LAYOUT_TABLE)
            if (!exists(mutationsFile))
            {
                log.error "TUMOUR_MUTATIONS_CSV file ${mutationsFile.toAbsolutePath()} does not exist."
                errors = true
            }
        }

        if (!containsKey('FASTA_REFERENCE'))
        {
            log.error "FASTA_REFERENCE must be set in your project's configuration file."
            errors = true
        }
        else
        {
            def refFile = file(FASTA_REFERENCE)
            if (!exists(refFile))
            {
                log.error "FASTA_REFERENCE file ${refFile.toAbsolutePath()} does not exist."
                errors = true
            }
        }

        if (!containsKey('THOUSAND_GENOMES_DATABASE'))
        {
            log.error "THOUSAND_GENOMES_DATABASE must be set in your project's configuration file."
            errors = true
        }
        else
        {
            def refFile = file(THOUSAND_GENOMES_DATABASE)
            if (!exists(refFile))
            {
                log.error "THOUSAND_GENOMES_DATABASE file ${refFile.toAbsolutePath()} does not exist."
                errors = true
            }
        }

        if (!containsKey('COSMIC_DATABASE'))
        {
            log.error "COSMIC_DATABASE must be set in your project's configuration file."
            errors = true
        }
        else
        {
            def refFile = file(COSMIC_DATABASE)
            if (!exists(refFile))
            {
                log.error "COSMIC_DATABASE file ${refFile.toAbsolutePath()} does not exist."
                errors = true
            }
        }

        /**
         * Setting values
         */

        if (isInteger(MAPPING_QUALITY))
        {
            if (MAPPING_QUALITY < 0)
            {
                log.error "MAPPING_QUALITY must be a positive integer."
                errors = true
            }
        }
        else
        {
            log.error "MAPPING_QUALITY must be an integer."
            errors = true
        }

        if (isInteger(BASE_QUALITY))
        {
            if (BASE_QUALITY < 1)
            {
                log.error "BASE_QUALITY must be a positive integer."
                errors = true
            }
        }
        else
        {
            log.error "BASE_QUALITY must be an integer."
            errors = true
        }

        if (isInteger(MPILEUP_MINIMUM_DEPTH))
        {
            if (MPILEUP_MINIMUM_DEPTH < 1)
            {
                log.error "MPILEUP_MINIMUM_DEPTH must be a positive integer > 0."
                errors = true
            }
        }
        else
        {
            log.error "MPILEUP_MINIMUM_DEPTH must be an integer."
            errors = true
        }

        if (isInteger(SLOP_BASES))
        {
            if (SLOP_BASES < 0)
            {
                log.error "SLOP_BASES must be a positive integer or zero."
                errors = true
            }
        }
        else
        {
            log.error "SLOP_BASES must be an integer."
            errors = true
        }

        if (!isBoolean(REMOVE_DUPLICATES))
        {
            log.error "REMOVE_DUPLICATES must be true or false."
            errors = true
        }

        if (isInteger(MAXIMUM_DEPTH))
        {
            if (MAXIMUM_DEPTH < 1)
            {
                log.error "MAXIMUM_DEPTH must be a positive integer > 0."
                errors = true
            }
        }
        else
        {
            log.error "MAXIMUM_DEPTH must be an integer."
            errors = true
        }

        if (isInteger(MINIMUM_REFERENCE_DEPTH))
        {
            if (MINIMUM_REFERENCE_DEPTH < 1)
            {
                log.error "MINIMUM_REFERENCE_DEPTH must be a positive integer > 0."
                errors = true
            }
        }
        else
        {
            log.error "MINIMUM_REFERENCE_DEPTH must be an integer."
            errors = true
        }

        if (isNumber(MQSB_THRESHOLD))
        {
            if (MQSB_THRESHOLD < 0d)
            {
                log.error "MQSB_THRESHOLD must be a positive number or zero."
                errors = true
            }
        }
        else
        {
            log.error "MQSB_THRESHOLD must be a number."
            errors = true
        }

        if (isInteger(ALT_ALLELES_THRESHOLD))
        {
            if (ALT_ALLELES_THRESHOLD < 0)
            {
                log.error "ALT_ALLELES_THRESHOLD must be a positive integer or zero."
                errors = true
            }
        }
        else
        {
            log.error "ALT_ALLELES_THRESHOLD must be an integer."
            errors = true
        }

        if (isInteger(MINOR_ALT_ALLELE_THRESHOLD))
        {
            if (MINOR_ALT_ALLELE_THRESHOLD < 0)
            {
                log.error "MINOR_ALT_ALLELE_THRESHOLD must be a positive integer or zero."
                errors = true
            }
        }
        else
        {
            log.error "MINOR_ALT_ALLELE_THRESHOLD must be an integer."
            errors = true
        }

        if (isInteger(COSMIC_THRESHOLD))
        {
            if (COSMIC_THRESHOLD < 0)
            {
                log.error "COSMIC_THRESHOLD must be a positive integer or zero."
                errors = true
            }
        }
        else
        {
            log.error "COSMIC_THRESHOLD must be an integer."
            errors = true
        }

        if (isNumber(PROPORTION_OF_CONTROLS))
        {
            if (PROPORTION_OF_CONTROLS < 0d)
            {
                log.error "PROPORTION_OF_CONTROLS must be a positive number or zero."
                errors = true
            }
        }
        else
        {
            log.error "PROPORTION_OF_CONTROLS must be a number."
            errors = true
        }

        if (isNumber(MAXIMUM_BACKGROUND_MEAN_ALLELE_FREQUENCY))
        {
            if (MAXIMUM_BACKGROUND_MEAN_ALLELE_FREQUENCY <= 0d)
            {
                log.error "MAXIMUM_BACKGROUND_MEAN_ALLELE_FREQUENCY must be a positive number."
                errors = true
            }
        }
        else
        {
            log.error "MAXIMUM_BACKGROUND_MEAN_ALLELE_FREQUENCY must be a number."
            errors = true
        }

        if (!isBoolean(IS_BLOODSPOT))
        {
            log.error "IS_BLOODSPOT must be true or false."
            errors = true
        }

        if (isNumber(OUTLIER_SUPPRESSION_THRESHOLD))
        {
            if (OUTLIER_SUPPRESSION_THRESHOLD < 0d)
            {
                log.error "OUTLIER_SUPPRESSION_THRESHOLD must be a positive number or zero."
                errors = true
            }
        }
        else
        {
            log.error "OUTLIER_SUPPRESSION_THRESHOLD must be a number."
            errors = true
        }

        if (isInteger(MINIMUM_FRAGMENT_LENGTH))
        {
            if (MINIMUM_FRAGMENT_LENGTH < 0)
            {
                log.error "MINIMUM_FRAGMENT_LENGTH must be a positive integer or zero."
                errors = true
            }
        }
        else
        {
            log.error "MINIMUM_FRAGMENT_LENGTH must be an integer."
            errors = true
        }

        if (isInteger(MAXIMUM_FRAGMENT_LENGTH))
        {
            if (MAXIMUM_FRAGMENT_LENGTH < 0)
            {
                log.error "MINIMUM_FRAGMENT_LENGTH must be a positive integer or zero."
                errors = true
            }
        }
        else
        {
            log.error "MAXIMUM_FRAGMENT_LENGTH must be an integer."
            errors = true
        }

        if (isInteger(MINIMUM_FRAGMENT_LENGTH) && isInteger(MAXIMUM_FRAGMENT_LENGTH) &&
            MINIMUM_FRAGMENT_LENGTH > MAXIMUM_FRAGMENT_LENGTH)
        {
            log.error "MINIMUM_FRAGMENT_LENGTH must be <= MAXIMUM_FRAGMENT_LENGTH"
            errors = true
        }

        if (isNumber(SMOOTHING))
        {
            if (SMOOTHING < 0d)
            {
                log.error "SMOOTHING must be a positive number or zero."
                errors = true
            }
        }
        else
        {
            log.error "SMOOTHING must be a number."
            errors = true
        }

        if (!isBoolean(ONLY_WEIGH_MUTANTS))
        {
            log.error "ONLY_WEIGH_MUTANTS must be true or false."
            errors = true
        }
    }

    return !errors
}

def validateReferenceFiles(params)
{
    def errors = false
    def layoutFile = file(params.LAYOUT_TABLE)
    def tumourMutationsFile = file(params.TUMOUR_MUTATIONS_CSV)

    log.info "Set to process files in the ${params.STUDY} study/project."

    /*
     * Check the layout file.
     */

    def requiredColumns = [ 'STUDY', 'SAMPLE_ID', 'BAM_FILE', 'CASE_OR_CONTROL', 'PATIENT',
                            'INPUT_INTO_LIBRARY_NG', 'SAMPLE_NAME', 'SAMPLE_TYPE', 'TIMEPOINT' ]
    def activeMarkers = [ '', 'yes', 'y', 'true', 't' ]

    layoutFile.withReader
    {
        reader ->

        def layoutContent = com.xlson.groovycsv.CsvParser.parseCsv(reader)
        def lineNumber = 1
        def uniqueSampleIds = new HashSet()
        def uniqueBamFiles = new HashSet()
        def activeBamCount = 0

        for (line in layoutContent)
        {
            if (lineNumber == 1)
            {
                // First time through, check columns.

                for (col in requiredColumns)
                {
                    if (!line.columns.containsKey(col))
                    {
                        log.error "${col} column is missing from the layout file."
                        errors = true
                    }
                }

                if (errors)
                {
                    break
                }
            }

            ++lineNumber

            if (isBlank(line.STUDY))
            {
                log.error "No STUDY on line ${lineNumber}."
                errors = true
            }

            if (isBlank(line.SAMPLE_ID))
            {
                log.error "No SAMPLE_ID on line ${lineNumber}."
                errors = true
            }
            else if (line.STUDY == params.STUDY)
            {
                if (uniqueSampleIds.contains(line.SAMPLE_ID))
                {
                    log.error "Have a duplicate SAMPLE_ID \"${line.SAMPLE_ID}\" within ${params.STUDY} (line ${lineNumber})."
                    errors = true
                }
                uniqueSampleIds << line.SAMPLE_ID
            }

            if (line.STUDY == params.STUDY)
            {
                def active = !line.columns.containsKey('ACTIVE') || trimToEmpty(line.ACTIVE).toLowerCase() in activeMarkers

                if (isBlank(line.BAM_FILE))
                {
                    log.error "There is no BAM_FILE for sample ${line.SAMPLE_ID} (line ${lineNumber})."
                    errors = true
                }
                else
                {
                    if (uniqueBamFiles.contains(line.BAM_FILE))
                    {
                        log.error "Have a duplicate BAM_FILE within ${params.STUDY} (line ${lineNumber})."
                        errors = true
                    }
                    else
                    {
                        if (active)
                        {
                            def foundBamFile = findFileOnPath(params.BAM_PATH, line.BAM_FILE)
                            if (foundBamFile)
                            {
                                ++activeBamCount
                            }
                            else
                            {
                                log.error "Could not locate \"${line.BAM_FILE}\" on BAM_PATH."
                                errors = true
                            }
                        }
                    }
                    uniqueBamFiles << line.BAM_FILE
                }

                if (active)
                {
                    try
                    {
                        Double.valueOf(line.INPUT_INTO_LIBRARY_NG)
                    }
                    catch (NumberFormatException e)
                    {
                        log.error "INPUT_INTO_LIBRARY_NG is not a number (sample ${line.SAMPLE_ID}, line ${lineNumber})."
                        errors = true
                    }
                }
            }
        }

        if (lineNumber == 1)
        {
            log.warn "${params.LAYOUT_TABLE} is empty."
            errors = true
        }

        if (activeBamCount == 0)
        {
            log.warn "There are no samples active for ${params.STUDY}."
            errors = true
        }
    }

    /*
     * Check the tumour mutations file.
     * We don't check the rows of this file, just that the required columns are present.
     */

    requiredColumns = [ 'CHROM', 'POS', 'REF', 'ALT', 'TUMOUR_AF', 'PATIENT' ]

    tumourMutationsFile.withReader
    {
        reader ->

        def tmIterator = com.xlson.groovycsv.CsvParser.parseCsv(reader)

        if (tmIterator.hasNext())
        {
            def line = tmIterator.next()

            for (col in requiredColumns)
            {
                if (!line.columns.containsKey(col))
                {
                    log.error "${col} column is missing from the tumour mutations file."
                    errors = true
                }
            }
        }
        else
        {
            log.warn "${params.TUMOUR_MUTATIONS_CSV} is empty."
            errors = true
        }
    }

    return !errors
}
