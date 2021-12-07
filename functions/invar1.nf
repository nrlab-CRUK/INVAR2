@Grab('org.apache.commons:commons-csv:1.9.0')

import org.apache.commons.csv.*

def checkBedFile(bedFile)
{
    file(bedFile).withReader
    {
        reader ->

        def message = "Your BED file is not of 1bp positions, aborting. Check your bed file\n" +
                      "NEEDS TO BE: CHR  POS-1 POS REF ALT" +
                      "No more than 1 line of header! Separated\tby\ttabs!"

        def line
        while (line = reader.readLine())
        {
            def parts = line.split("\t")
            if (parts.length != 5)
            {
                throw new Exception(message)
            }
            def start = parts[1] as long
            def end = parts[2] as long

            if (end - start != -1L)
            {
                throw new Exception(message)
            }
        }
    }
}
