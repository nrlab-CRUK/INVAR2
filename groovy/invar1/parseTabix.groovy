@Grab('org.apache.commons:commons-lang3:3.9')

import groovy.json.JsonSlurper
import static org.apache.commons.lang3.StringUtils.trimToEmpty

assert args.length == 5 : "Expect exactly five arguments for parseTabix.groovy"

final def mutationJson = args[0]
final def snpFile = new File(args[1])
final def cosmicFile = new File(args[2])
final def trinucleotideFile = new File(args[3])
final def outFile = new File(args[4])

def vcfParts =
{
    vcfLine ->
    
    def mainParts = vcfLine.split('\t')
    def subparts = mainParts[7].split(';')
    def info = [:]
    for (subpart in subparts)
    {
        def kv = subpart.split('=')
        def k = kv[0]
        def v = kv.length > 1 ? kv[1] : true
        info[k] = v
    }
    return info
}

def mutationInfo = new JsonSlurper().parseText(mutationJson)

def k1g_af = 0d
snpFile.withReader
{
    reader ->
    def vcfLine = reader.readLine()
    if (vcfLine)
    {
        def info = vcfParts(vcfLine)
        k1g_af = info['AF'].split(',')[0] as double
    }
}

cosmic_sampleN = 0
cosmic_snp = 0

cosmicFile.withReader
{
    reader ->
    def vcfLine = reader.readLine()
    if (vcfLine)
    {
        def info = vcfParts(vcfLine)
        cosmic_sampleN = info['CNT'] as int
        cosmic_snp = info.containsKey('SNP') ? 1 : 0
    }
}

trinucleotide = ""

trinucleotideFile.withReader
{
    reader ->
    trinucleotide = trimToEmpty(reader.readLine())
}

outFile.withPrintWriter
{
    pw ->
    pw.println("${mutationInfo.keySet().join('\t')}\tCOSMIC_MUTATIONS\tCOSMIC_SNP\t1KG_AF\tTRINUCLEOTIDE")
    pw.println("${mutationInfo.values().join('\t')}\t${cosmic_sampleN}\t${cosmic_snp}\t${k1g_af}\t${trinucleotide}")
}
