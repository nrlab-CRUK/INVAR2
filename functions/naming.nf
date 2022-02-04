/**
 * Convert a string into a safe version for a file name. Spaces are replaced
 * by a single underscore then any character that is not a word character
 * (letter, digit or underscore) is removed.
 */
def makeSafeForFileName(str)
{
    str.replaceAll(/\s+/, '_').replace(/[^\w]+/, '')
}
