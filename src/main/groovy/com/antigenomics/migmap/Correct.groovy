/*
 * Copyright 2014-2015 Mikhail Shugay
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.migmap

import com.milaboratory.core.sequence.NucleotideSequence
import com.milaboratory.core.tree.SequenceTreeMap
import com.milaboratory.core.tree.TreeSearchParameters
import groovy.transform.Canonical

def DEFAULT_ERROR_RATE = "0.01", DEFAULT_DEPTH = "2"

def cli = new CliBuilder(usage: "Correct [options] clonotype_table.txt output.txt")

cli.r(args: 1, argName: "float", "Child-to-parent ratio for a single error. [default=$DEFAULT_ERROR_RATE]")
cli.d(args: 1, argName: "int", "Scan depth, maximum number of allowed mismatches. [default=$DEFAULT_DEPTH]")

// Misc
cli.h("Display this help message")

// PARSE ARGUMENTS

def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() != 2) {
    cli.usage()
    System.exit(3)
}

// I/O

def inputFileName = opt.arguments()[0], outputFileName = opt.arguments()[1]

if (!new File(inputFileName).exists()) {
    Util.error("Input file $inputFileName does not exist.", 3)
}

// Ensure all subdirs are created for output
def parentFolder = new File(outputFileName).absoluteFile.parentFile
if (parentFolder) {
    parentFolder.mkdirs()
}

def errorRate = Math.log((opt.r ?: DEFAULT_ERROR_RATE).toFloat()), depth = (opt.d ?: DEFAULT_DEPTH).toInteger()

// Parsing required columns, here we operate on full contigs

def countColIndex = -1, freqColIndex = -1, contigNtIndex = -1

def header = (String[]) null
def parseColumns = {
    freqColIndex = header.findIndexOf { it.toLowerCase() == "freq" }
    countColIndex = header.findIndexOf { it.toLowerCase() == "count" }
    contigNtIndex = header.findIndexOf { it.toLowerCase() == "contignt" }

    if ([freqColIndex, countColIndex, contigNtIndex].any { it < 0 }) {
        Util.error("One of the critical columns ('freq', 'count', 'contignt') is missing.", 3)
    }
}

@Canonical
class ClonotypeEntry {
    int count
    float freq
    String[] data
    NucleotideSequence seq
    ClonotypeEntry parent

    void append(ClonotypeEntry other) {
        count += other.count
        freq += other.freq
    }
}

def stm = new SequenceTreeMap<NucleotideSequence, ClonotypeEntry>(NucleotideSequence.ALPHABET)
def searchParams = new TreeSearchParameters(depth, depth, depth, depth)

new File(inputFileName).splitEachLine("\t") { String[] splitLine ->
    if (!header) {
        header = splitLine
        parseColumns()
    } else {
        def seq = new NucleotideSequence(splitLine[contigNtIndex])
        stm.put(seq, new ClonotypeEntry(splitLine[countColIndex].toInteger(),
                splitLine[freqColIndex].toFloat(),
                splitLine, seq))
    }
}

stm.values().each { child ->
    def iter = stm.getNeighborhoodIterator(child.seq, searchParams)
    float bestParentScore = 1.0

    def bestParent = null, parent

    while ((parent = iter.next())) {
        if (parent.count > child.count) {
            def score = Math.log(child.count / (float) parent.count) / iter.mutationsCount

            if (score < bestParentScore) {
                bestParent = parent
                bestParentScore = score
            }
        }
    }

    if (bestParentScore <= errorRate) {
        child.parent = bestParent
    }
}

stm.values().sort { it.count }.each {
    if (it.parent != null) {
        it.parent.append(it)
    }
}

new File(outputFileName).withPrintWriter { pw ->
    pw.println(header.join("\t"))

    stm.values().each {
        if (it.parent == null) {
            pw.println(it.data.join("\t"))
        }
    }
}


