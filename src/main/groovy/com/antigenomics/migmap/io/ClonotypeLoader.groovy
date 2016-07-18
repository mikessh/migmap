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

package com.antigenomics.migmap.io

import com.antigenomics.migmap.blast.PSegments
import com.antigenomics.migmap.clonotype.Clonotype
import com.antigenomics.migmap.genomic.Segment
import com.antigenomics.migmap.genomic.SegmentType
import com.antigenomics.migmap.mapping.Cdr3Markup
import com.antigenomics.migmap.mapping.Truncations
import com.antigenomics.migmap.mutation.Mutation
import com.antigenomics.migmap.mutation.SubRegion
import groovy.transform.CompileStatic

@CompileStatic
class ClonotypeLoader {
    static final int cdr3ntInd, cdr3aaInd,
                     vInd, dInd, jInd,
                     countInd, freqInd,
                     hasCdr3Ind, inFrameInd, noStopInd, completeInd, canonicalInd,
                     vEndInCdr3Ind, dStartInCdr3Ind, dEndInCdr3Ind, jStartInCdr3Ind,
                     vDelInd, dDel5Ind, dDel3Ind, jDelInd,
                     polVInd, polD5Ind, polD3Ind, polJInd
    static final Map<String, Integer> columnIndexMap = new HashMap<String, Integer>()

    static {
        Clonotype.OUTPUT_HEADER.split("\t").eachWithIndex { String entry, int i ->
            columnIndexMap.put(entry, i)
        }

        cdr3ntInd = columnIndexMap["cdr3nt"]
        cdr3aaInd = columnIndexMap["cdr3aa"]
        vInd = columnIndexMap["v"]
        dInd = columnIndexMap["d"]
        jInd = columnIndexMap["j"]
        countInd = columnIndexMap["count"]
        freqInd = columnIndexMap["freq"]
        hasCdr3Ind = columnIndexMap["has.cdr3"]
        inFrameInd = columnIndexMap["in.frame"]
        noStopInd = columnIndexMap["no.stop"]
        completeInd = columnIndexMap["complete"]
        canonicalInd = columnIndexMap["canonical"]
        vEndInCdr3Ind = columnIndexMap["v.end.in.cdr3"]
        dStartInCdr3Ind = columnIndexMap["d.start.in.cdr3"]
        dEndInCdr3Ind = columnIndexMap["d.end.in.cdr3"]
        jStartInCdr3Ind = columnIndexMap["j.start.in.cdr3"]
        vDelInd = columnIndexMap["v.del"]
        dDel5Ind = columnIndexMap["d.del.3"]
        dDel3Ind = columnIndexMap["d.del.5"]
        jDelInd = columnIndexMap["j.del"]
        polVInd = columnIndexMap["pol.v"]
        polD5Ind = columnIndexMap["pol.d.5"]
        polD3Ind = columnIndexMap["pol.d.3"]
        polJInd = columnIndexMap["pol.j"]
    }

    static List<Clonotype> load(String plainTextOutput) {
        def clonotypes = new ArrayList<Clonotype>()
        boolean firstLine = true

        new File(plainTextOutput).splitEachLine("\t") { splitLine ->
            if (firstLine) {
                if (!splitLine.join("\t").startsWith(Clonotype.OUTPUT_HEADER)) {
                    throw new RuntimeException("Bad clonotype table header")
                }
                firstLine = false
            } else {
                def clonotype = new Clonotype(
                        splitLine[cdr3ntInd],
                        splitLine[cdr3aaInd],
                        decodeSegment(splitLine[vInd]),
                        decodeSegment(splitLine[dInd]),
                        decodeSegment(splitLine[jInd]),
                        decodeMutations(splitLine),
                        splitLine[countInd].toLong(),
                        splitLine[freqInd].toDouble(),
                        null,
                        null,
                        splitLine[hasCdr3Ind].toBoolean(),
                        splitLine[inFrameInd].toBoolean(),
                        splitLine[noStopInd].toBoolean(),
                        splitLine[completeInd].toBoolean(),
                        splitLine[canonicalInd].toBoolean(),
                        decodeCdr3Markup(splitLine),
                        decodeTruncations(splitLine),
                        decodePSegments(splitLine)
                )

                clonotypes.add(clonotype)
            }
        }

        clonotypes
    }

    static List<Mutation> decodeMutations(String[] splitLine) {
        def mutations = new ArrayList<Mutation>()
        SubRegion.REGION_LIST.each { SubRegion subRegion ->
            def ntMutations = splitLine[columnIndexMap["mutations.nt." + subRegion.toString()]].split(","),
                aaMutations = splitLine[columnIndexMap["mutations.aa." + subRegion.toString()]].split(",")

            for (int i = 0; i < ntMutations.length; i++) {
                def ntMutation = ntMutations[i],
                    aaMutation = aaMutations[i]

                def mutation = Mutation.fromString(ntMutation, aaMutation)

                mutation.subRegion = subRegion
                // todo: Parent segment (need segment database and some special coding magick)
            }
        }
        mutations
    }

    static Cdr3Markup decodeCdr3Markup(String[] splitLine) {
        new Cdr3Markup(splitLine[vEndInCdr3Ind].toInteger(),
                splitLine[dStartInCdr3Ind].toInteger(),
                splitLine[dEndInCdr3Ind].toInteger(),
                splitLine[jStartInCdr3Ind].toInteger())
    }

    static Truncations decodeTruncations(String[] splitLine) {
        new Truncations(splitLine[vDelInd].toInteger(),
                splitLine[dDel5Ind].toInteger(),
                splitLine[dDel3Ind].toInteger(),
                splitLine[jDelInd].toInteger())
    }

    static PSegments decodePSegments(String[] splitLine) {
        new PSegments(splitLine[polVInd].toInteger(),
                splitLine[polD5Ind].toInteger(),
                splitLine[polD3Ind].toInteger(),
                splitLine[polJInd].toInteger())
    }

    static Segment decodeSegment(String name) {
        // todo: use SegmentDatabase
        new Segment(null, SegmentType.valueOf(name[3]), name[0..2], name, null, -1)
    }
}
