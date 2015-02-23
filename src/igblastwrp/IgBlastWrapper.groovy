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

def cli = new CliBuilder(usage: 'igblastwrp [options] input.(fa/fastq)[.gz] outputPrefix')
cli.h('usage')

cli.f('Report functional (no stop codons and frameshifts) clonotypes only')
cli.c('Report clonotypes with complete CDR3s only')
cli._(longOpt: 'no-cdr3', 'Report clonotypes for which no CDR3 was identified')

cli.l(args: 1, argName: '0, 1 or 2',
        'Level of clonotype detalization: 0 - CDR3, 1 - CDR1,2,3, 2 - CDR1,2,3+V-mutations. ' +
                'Several levels could be specified as comma-separated string [default = 2]. ' +
                'Also specifies stop/functional filtering scope.')

cli.R(args: 1, argName: 'TRA|B|G|D and IGH|L|K', 'Receptor gene and chain, e.g. \'TRA\' [required]')
cli._(longOpt: 'all-alleles', 'Use all alleles during alignment (this is going to be slower). ' +
        'Will use only major (*01) alleles if option is not set.')
cli._(longOpt: 'dump-unmapped', args: 1, argName: 'path',
        'Specifies a file that output sequences that were not mapped to V(D)J segments.')
cli.a('[migec-compatibility] Assume MIG-assembled data. Header should contain UMI:NNNNN:count')
cli.S(args: 1, argName: 'species name', 'human, mouse, rat, rabbit or rhesus_monkey [default=human]')
cli.q(args: 1, 'quality threshold, 2..40 [default = 25]')
cli.p(args: 1, 'number of threads to use [default = all available processors]')
cli._(longOpt: 'data-dir', args: 1, argName: 'path', 'path to folder that contains IgBlastWrapper bundle: ' +
        'data/ and bin/ folders [default = parent directory of script]')
cli._(longOpt: "debug", "Debug mode. Will run in a single thread and pring IgBlast output to stdout.")
cli.h('display help message')

cli.N(args: 1, 'number of reads to take')

def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() < 2 || !opt.R) {
    cli.usage()
    System.exit(0)
}

def levels = (opt.l ?: "2").split(",").collect { it.toInteger() }.unique()
levels.each { level ->
    if (level < 0 || level > 2) {
        println "ERROR Illegal clonotype detalization level $level"
        System.exit(-1)
    }
}
int N = (opt.N ?: "-1").toInteger()
String SPECIES = opt.S ?: "human", GENE = opt.R[0..1], CHAIN = opt.R[2]
String outUnmapped = opt.'dump-unmapped'
byte qualThreshold = (opt.q ?: "25").toInteger()
boolean funcOnly = opt.f, completeOnly = opt.c, assembled = opt.a, reportNoCdr3 = opt.'no-cdr3',
        allAlleles = opt.'all-alleles', debug = opt.'debug'
int THREADS = debug ? 1 : ((opt.p ?: "${Runtime.runtime.availableProcessors()}").toInteger())

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

String inputFileName = opt.arguments()[0], outputFilePrefix = opt.arguments()[1]

if (!new File(inputFileName).exists()) {
    println "[ERROR] Input file doesn't exist"
    System.exit(-1)
}

def outputDir = new File(outputFilePrefix).absoluteFile.parentFile
if (!outputDir.exists()) {
    println "[INFO] Output path doesn't exist.. creating"
    def created = outputDir.mkdirs()
    if (!created) {
        println "[ERROR] Failed to create output path"
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

println "[${new Date()}] Reading $inputFileName.."
Read read
int n = 0
while ((read = reader.next()) != null) {
    def seqData = nonRedundantSequenceMap.get(read.seq)

    if (!seqData)
        nonRedundantSequenceMap.put(read.seq, seqData = new SeqData(nonRedundantSequenceMap.size(), read.qual))

    int count = 1
    if (assembled) {
        def umiFieldMatcher = read.header =~ /UMI:(?:[ATGC]+):([0-9]+)/
        if (umiFieldMatcher.size() > 0)
            count = (umiFieldMatcher[0][1].toString()).toInteger()
        else {
            println "[ERROR] Read header misses UMI field in assembly mode"
            System.exit(-1)
        }
    }
    seqData.append(read.qual, count)

    n++

    if (n % 50000 == 0)
        println "[${new Date()}] $n sequences read, ${nonRedundantSequenceMap.size()} non-redundant ones so far.."

    if (N > 0 && n == N)
        break
}
println "[${new Date()}] Finished reading, $n sequences total, ${nonRedundantSequenceMap.size()} non-redundant ones"
println "[${new Date()}] Preparing to run IgBlast (creating .fa chunks)"

//
// Create .fa chunks for IgBlast
//
def fastaChunks = new ArrayList<String>()
def seqIter = nonRedundantSequenceMap.entrySet().iterator()
for (int i = 0; i < THREADS; i++) {
    def prefix = "_igblastwrp-" + UUID.randomUUID().toString() + "_" + i
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

println "[${new Date()}] Running IgBlast for $inputFileName and parsing output"

def listener = Thread.start {
    while (!finished) {
        Thread.sleep(30000)

        int m = clonotypeMap.size(), M = nonRedundantSequenceMap.size()
        println "[${new Date()}] $m non-redundant sequences succesfully " +
                "aligned and processed so far " +
                "(${m > 0 ? ((int) (10000 * (double) m / (double) M)) / 100 : 0}%)"
    }
}

if (debug) {
    new BlastRunner(THREADS, SCRIPT_PATH, SPECIES, GENE, CHAIN, allAlleles, fastaChunks[0], clonotypeMap).runIdle()
    System.exit(0)
} else {
    if (THREADS > 1) {
        def pool = Executors.newFixedThreadPool(THREADS)

        (0..(THREADS - 1)).collect { p ->
            pool.submit(new BlastRunner(THREADS, SCRIPT_PATH, SPECIES, GENE, CHAIN, allAlleles, fastaChunks[p], clonotypeMap))
        }
        pool.shutdown()

        while (!pool.terminated);
    } else {
        new BlastRunner(THREADS, SCRIPT_PATH, SPECIES, GENE, CHAIN, allAlleles, fastaChunks[0], clonotypeMap).run()
    }
}

listener.setDefaultUncaughtExceptionHandler({ t, ex ->
    int m = clonotypeMap.size(), M = nonRedundantSequenceMap.size()
    println "[${new Date()}] Finished. $m non-redundant sequences succesfully " +
            "aligned and processed " +
            "(${m > 0 ? ((int) (10000 * (double) m / (double) M)) / 100 : 0}%)"
} as Thread.UncaughtExceptionHandler)
finished = true
listener.interrupt()
listener.join()

//
// Group clonotypes
//
def dumpUnmapped = outUnmapped != null
def pwDump = dumpUnmapped ? new File(outUnmapped).newPrintWriter() : null
levels.each { level ->
    def resultsMap = new HashMap<String, ClonotypeData>()
    def outputFileName = outputFilePrefix + ".L${level}.txt"

    println "[${new Date()}] Generating clonotypes for level $level"
    nonRedundantSequenceMap.eachWithIndex { it, ind ->
        def clonotype = clonotypeMap[it.value.seqId.toString()]
        if (clonotype != null) {
            def clonotypeKey = clonotype.generateKey(it.key, level, funcOnly, completeOnly, reportNoCdr3)
            if (clonotypeKey != null) {
                def clonotypeData = resultsMap[clonotypeKey]

                def qual = it.value.computeQual()

                if (!clonotypeData)
                    resultsMap.put(clonotypeKey, clonotype.appendToData(null, qual, qualThreshold,
                            it.value.nReads, it.value.nEvents, level))  // create new
                else
                    clonotype.appendToData(clonotypeData, qual, qualThreshold,
                            it.value.nReads, it.value.nEvents, level)
            } else if (dumpUnmapped) {
                pwDump.println("INCOMPLETE SEQ$ind READS:${it.value.nReads} EVENTS:${it.value.nEvents}\n" +
                        it.key +
                        "\n+\n" +
                        it.value.computeQual())
            }
        } else if (dumpUnmapped) {
            pwDump.println("UNMAPPED SEQ$ind READS:${it.value.nReads} EVENTS:${it.value.nEvents}\n" +
                    it.key +
                    "\n+\n" +
                    it.value.computeQual())
        }
    }

    if (pwDump)
        pwDump.close()

//
// Write output
//
    println "[${new Date()}] Writing output to $outputFileName"
    def nReads = 0, nEvents = 0
    new File(outputFileName).withPrintWriter { pw ->
        pw.println "#reads_count\treads_freq\tmig_count\tmig_freq\t" +
                Clonotype.KEY_HEADER + "\t" + ClonotypeData.VALUE_HEADER

        // quality filter (CDR3)
        resultsMap = resultsMap.findAll {
            byte minQual = it.value.summarizeQuality()
            minQual >= qualThreshold
        }

        // count totals
        resultsMap.each {
            nReads += it.value.nReads
            nEvents += it.value.nEvents
        }

        // store clonotype table
        resultsMap.sort { -it.value.nReads }.each {
            def key = it.key
            key = key.substring(0, key.lastIndexOf("\t")) // trim hypermutations, they will be taken from 'value'

            pw.println(it.value.nReads + "\t" + (it.value.nReads / nReads) + "\t" +
                    it.value.nEvents + "\t" + (it.value.nEvents / nEvents) + "\t" +
                    key + "\t" + it.value.toString())
        }
    }
}

println "[${new Date()}] Finished"