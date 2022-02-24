/**
 * Validation pipeline configuration.
 * Make sure mandatory parameters are defined and all others are within acceptable ranges.
 */

import static java.nio.file.Files.exists

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

def validatePipeline(params)
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

        def inputListFile = file(INPUT_FILES)
        if (!exists(inputListFile))
        {
            log.error "INPUT_FILES file ${inputListFile.toAbsolutePath()} does not exist."
            errors = true
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

        if (isInteger(MPILEUP_MAXIMUM_DEPTH))
        {
            if (MPILEUP_MAXIMUM_DEPTH < 1)
            {
                log.error "MPILEUP_MAXIMUM_DEPTH must be a positive integer > 0."
                errors = true
            }
        }
        else
        {
            log.error "MPILEUP_MAXIMUM_DEPTH must be an integer."
            errors = true
        }

        if (isInteger(MPILEUP_MINIMUM_DEPTH) && isInteger(MPILEUP_MAXIMUM_DEPTH) &&
            MPILEUP_MINIMUM_DEPTH > MPILEUP_MAXIMUM_DEPTH)
        {
            log.error "MPILEUP_MINIMUM_DEPTH must be <= MPILEUP_MAXIMUM_DEPTH"
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
