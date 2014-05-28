package igblastwrp

import igblastwrp.blast.BlastRunner
import igblastwrp.blast.Clonotype
import igblastwrp.blast.ClonotypeData
import igblastwrp.io.FastaReader
import igblastwrp.io.FastqReader
import igblastwrp.io.Read
import igblastwrp.io.SeqData

import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.Executors

/**
 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

def cli = new CliBuilder(usage: 'igblastwrp [options] input.(fa/fastq)[.gz] output')
cli.h('usage')

cli.f('Report clonotypes with functional CDR3s only')
cli.c('Report clonotypes with complete CDR3s only')

cli.R(args: 1, argName: 'TRA|B|G|D and IGH|L|K', 'Receptor gene and chain, e.g. \'TRA\' [required]')
cli.S(args: 1, argName: '\'human\' or \'mouse\'', 'Species [default=human]')
cli.q(args: 1, 'quality threshold, 2..40 [default = 25]')
cli.p(args: 1, 'number of threads to use [default = all available processors]')
cli._(longOpt: 'data-dir', args: 1, argName: 'path', 'path to folder that contains IgBlastWrapper bundle: ' +
        'data/ and bin/ folders [default = parent directory of script]')
cli.N(args: 1, 'number of reads to take')

def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() < 2 || !opt.R) {
    cli.usage()
    System.exit(0)
}

int THREADS = (opt.p ?: "${Runtime.runtime.availableProcessors()}").toInteger()
int N = (opt.N ?: "-1").toInteger()
String SPECIES = opt.S ?: "human", GENE = opt.R[0..1], CHAIN = opt.R[2]
byte qualThreshold = (opt.q ?: "25").toInteger()
boolean funcOnly = opt.f, completeOnly = opt.c

//
// DATA BUNDLE
//

String SCRIPT_PATH = opt.'data-dir' ? new File(opt.'data-dir'.toString()) : null

if (!SCRIPT_PATH) {
    def SCRIPT_SOURCE = new File(getClass().protectionDomain.codeSource.location.path)
    SCRIPT_PATH = SCRIPT_SOURCE.parent.replaceAll("%20", " ")

    if (SCRIPT_SOURCE.absolutePath.endsWith(".groovy")) // trim /src for script
        SCRIPT_PATH = SCRIPT_PATH.replaceAll(/(?:src\/){1}.+/, "")
} else {
    def scriptParentDir = new File(SCRIPT_PATH)
    if (!scriptParentDir.exists()) {
        println "Bad path to data bundle"
        System.exit(-1)
    }
    SCRIPT_PATH = scriptParentDir.absolutePath
}

//
// IO
//

String inputFileName = opt.arguments()[0], outputFileName = opt.arguments()[1]

if (!new File(inputFileName).exists()) {
    println "Input file doesn't exist"
    System.exit(-1)
}

def outputDir = new File(outputFileName).absoluteFile.parentFile
if (!outputDir.exists()) {
    println "Output path doesn't exist.. creating"
    def created = outputDir.mkdirs()
    if (!created) {
        println "Failed to create output path"
        System.exit(0)
    }
}

//
//
// SCRIPT BODY
//
//

//
// Read input, group reads
//
def nonRedundantSequenceMap = new HashMap<String, SeqData>()
def reader = inputFileName =~ /fastq(?:\.gz)?$/ ? new FastqReader(inputFileName) :
        new FastaReader(inputFileName)

println "[${new Date()}] Reading input.."
Read read
int n = 0
while ((read = reader.next()) != null) {
    def seqData = nonRedundantSequenceMap.get(read.seq)

    if (!seqData)
        nonRedundantSequenceMap.put(read.seq, seqData = new SeqData(nonRedundantSequenceMap.size(), read.qual))

    seqData.append(read.qual)

    n++

    if (n % 50000 == 0)
        println "[${new Date()}] $n sequences read, ${nonRedundantSequenceMap.size()} non-redundant ones so far.."

    if (N > 0 && n == N)
        break
}
println "[${new Date()}] Finished reading, $n sequences total, ${nonRedundantSequenceMap.size()} non-redundant ones"
println "[${new Date()}] Preparing to run IgBlast"

//
// Create .fa chunks for IgBlast
//
def fastaChunks = new ArrayList<String>()
def seqIter = nonRedundantSequenceMap.entrySet().iterator()
for (int i = 0; i < THREADS; i++) {
    def prefix = "_igblastwrp-" + UUID.randomUUID().toString()
    def chunkFileName = "$outputDir/${prefix}.fa"

    new File(chunkFileName).withPrintWriter { pw ->
        for (int j = 0; j < nonRedundantSequenceMap.size() / THREADS; j++) {
            if (!seqIter.hasNext())
                break

            def seqEntry = seqIter.next()

            pw.println(">$seqEntry.value.seqId")
            pw.println(seqEntry.key)
        }
    }
    new File(chunkFileName).deleteOnExit()
    fastaChunks.add(chunkFileName)
}

//
// Run IgBlast in parallel
//
def clonotypeMap = new ConcurrentHashMap<String, Clonotype>()
boolean finished = false

println "[${new Date()}] Running IgBlast and parsing output"

def listener = Thread.start {
    while (!finished) {
        sleep 30000
        println "[${new Date()}] ${clonotypeMap.size()} non-redundant sequences succesfully " +
                "aligned and processed so far " +
                "(${((int) (10000 * clonotypeMap.size() / nonRedundantSequenceMap.size())) / 100}%)"
    }
}

def pool = Executors.newFixedThreadPool(THREADS)
(0..(THREADS - 1)).collect { p ->
    pool.submit(new BlastRunner(THREADS, SCRIPT_PATH, SPECIES, GENE, CHAIN, fastaChunks[p], clonotypeMap))
}

pool.shutdown()

while (!pool.terminated);

finished = true
listener.join()

println "[${new Date()}] Finished. ${clonotypeMap.size()} non-redundant sequences succesfully " +
        "aligned and processed (${((int) (10000 * clonotypeMap.size() / nonRedundantSequenceMap.size())) / 100}%)"

//
// Group clonotypes
//
def resultsMap = new HashMap<String, ClonotypeData>()

println "[${new Date()}] Generating clonotypes"
nonRedundantSequenceMap.each {
    def clonotype = clonotypeMap[it.value.seqId.toString()]
    if (clonotype != null &&
            (!funcOnly || clonotype.functional) &&
            (!completeOnly || clonotype.complete)) {
        def clonotypeEntry = clonotype.generateEntry(it.key)
        def clonotypeData = resultsMap[clonotypeEntry]

        String qual = it.value.computeQual()

        if (!clonotypeData)
            resultsMap.put(clonotypeEntry, clonotype.appendToData(null, qual))  // create new
        else
            clonotype.appendToData(clonotypeData, qual)
    }
}

//
// Write output
//
println "[${new Date()}] Writing output"
new File(outputFileName).withPrintWriter { pw ->
    pw.println "Count\t" + Clonotype.HEADER + "\t" + ClonotypeData.HEADER
    resultsMap.sort { -it.value.count }.each {
        byte minQual = it.value.summarizeQuality()
        if (minQual >= qualThreshold)
            pw.println(it.value.count + "\t" + it.key + "\t" + it.value.qualString())
    }
}

println "[${new Date()}] Finished"