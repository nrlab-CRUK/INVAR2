import java.util.regex.Pattern

/**
 * Find a file on a Unix/Windows style path list. The path is
 * searched in order, and the first file found with an exact name
 * match is returned. If there is no file with the name on any
 * path, a dummy path that won't exist is returned ("/dev/null/<filename>").
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
            return file.absolutePath
        }
    }

    // Return a file path that can (probably) never exist.
    return "/dev/null/${filename}"
}
