package igblastwrp

import igblastwrp.blast.BlastRunner
import igblastwrp.blast.Clonotype
import igblastwrp.io.FastaReader
import igblastwrp.io.FastqReader
import igblastwrp.io.Read
import igblastwrp.io.SeqData

import java.util.concurrent.ConcurrentHashMap

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
cli.C(args: 1, argName: '\'TRA\', \'TRB\', \'TRG\', \'TRD\',  \'IGL\', \'IGK\' or \'IGH\'', 'Receptor chain [required]')
cli.S(args: 1, argName: '\'human\' or \'mouse\'', 'Species [default=HomoSapiens]')
cli.p(args: 1, 'number of threads to use [default = all available processors]')
cli.N(args: 1, 'number of reads to take')

def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() < 2 || !opt.C) {
    cli.usage()
    System.exit(0)
}

int THREADS = (opt.p ?: "${Runtime.runtime.availableProcessors()}").toInteger()
String SPECIES = opt.S ?: "human", GENE = opt.C[0..1], CHAIN = opt.C[2]
String inputFileName = opt.arguments()[0], outputFileName = opt.arguments()[1]

if (!new File(inputFileName).exists()) {
    println "Input file doesn't exist"
    System.exit(-1)
}

def outputDir = new File(outputFileName).parentFile
if (!outputDir) {
    println "Output path doesn't exist.. creating"
    def created = outputDir.mkdirs()
    if (!created) {
        println "Failed to create output path"
        System.exit(0)
    }
}

//
// Read input, group reads
//
def seqRedundMap = new HashMap<String, SeqData>()
def reader = inputFileName =~ /fastq(?:\.gz)?$/ ? new FastqReader(inputFileName) :
        new FastaReader(inputFileName)

Read read
int readId = 0
while ((read = reader.next()) != null) {
    def seqData = seqRedundMap.get(read.seq)

    if (!seqData)
        seqRedundMap.put(read.seq, seqData = new SeqData(seqRedundMap.size()))

    seqData.readIds.add(readId++)
}

//
// Create .fa chunks for IgBlast
//
def fastaChunks = new ArrayList<String>()
for (int i = 0; i < THREADS; i++) {
    def prefix = UUID.randomUUID().toString()
    def chunkFileName = "$outputDir/${prefix}.fa"
    new File(chunkFileName).withPrintWriter { pw ->
        seqRedundMap.each {
            pw.println(">$it.value.seqId")
            pw.println(it.key)
        }
    }
    fastaChunks.add(chunkFileName)
}

//
// Run IgBlast in parallel
//
def clonotypeMap = new ConcurrentHashMap<String, Clonotype>()

def processes = (0..(THREADS - 1)).collect { p ->
    new Thread(new BlastRunner(SPECIES, GENE, CHAIN, fastaChunks[p], clonotypeMap))
}

processes.each { it.run() }
processes.each { it.join() }

//runner.run()

println "SeqId\t" + Clonotype.HEADER_RAW
clonotypeMap.each {
    println it.key + "\t" + it.value
}

//SPECIES = "human", GENE = "TR", CHAIN = "A",
//INPUT = "/Users/mikesh/Programming/igblastwrp/test_tra.fa",
//OUTPUT = "/Users/mikesh/Programming/igblastwrp/results_a.txt"
