/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.higblast.blast

import com.antigenomics.higblast.genomic.Segment
import com.antigenomics.higblast.genomic.SegmentDatabase
import com.antigenomics.higblast.mapping.Cdr3Markup
import com.antigenomics.higblast.mapping.Mapping
import com.antigenomics.higblast.mapping.RegionMarkup
import com.antigenomics.higblast.mapping.Truncations
import com.antigenomics.higblast.mutation.MutationExtractor

import java.util.concurrent.atomic.AtomicInteger

import static com.antigenomics.higblast.Util.BLAST_NA
import static com.antigenomics.higblast.Util.groomMatch

class BlastParser {
    final AtomicInteger total = new AtomicInteger(),
                        noMatch = new AtomicInteger(),
                        vNotFound = new AtomicInteger(),
                        noCdr3 = new AtomicInteger()
    final SegmentDatabase segmentDatabase
    final JRefSearcher jRefSearcher

    public BlastParser(SegmentDatabase segmentDatabase) {
        this.segmentDatabase = segmentDatabase
        this.jRefSearcher = new JRefSearcher()
    }

    Mapping parse(String chunk) {
        total.incrementAndGet()

        def summary, alignments, cdrBounds

        // Rearrangement summary
        summary = segmentDatabase.hasD ?
                groomMatch(chunk =~
                        //         V     D     J  chain stop frame (prod) strand
                        /# V-.+\n(.+)\t(.+)\t(.+)\tV.\t(.+)\t(.+)\t.+\t(.+)/) :

                groomMatch(chunk =~
                        //         V     J    chain   stop frame (prod)  strand
                        /# V-.+\n(.+)\t(.+)\tV.\t(.+)\t(.+)\t.+\t(.+)/)

        if (summary == null) {
            noMatch.incrementAndGet()
            return
        }

        def rc = summary[-1] != "+", inFrame = summary[-2] != "Out-of-frame", noStop = summary[-3] != "Yes"

        // Information on segments mapped
        // - Segment names, can be multiple of them
        def vSegmentNames = summary[0],
            jSegmentNames = segmentDatabase.hasD ? summary[2] : summary[1],
            dSegmentNames = segmentDatabase.hasD ? summary[1] : BLAST_NA

        def dFound = dSegmentNames != BLAST_NA,
            vFound = vSegmentNames != BLAST_NA,
            jFound = jSegmentNames != BLAST_NA

        if (!vFound) {
            vNotFound.incrementAndGet()
            return null
        }

        def vSegment = segmentDatabase.segments[vSegmentNames.split(",")[0]],
            dSegment = dFound ? segmentDatabase.segments[dSegmentNames.split(",")[0]] : Segment.DUMMY_D,
            jSegment = jFound ? segmentDatabase.segments[jSegmentNames.split(",")[0]] : Segment.DUMMY_J

        // - Alignments for V, D and J segments, remember here and further BLAST coordinates are 1-based
        alignments = [
                groomMatch(chunk =~
                        //                          qstart     qseq        sstart     sseq
                        /# Hit table(?:.+\n)+V\t.+\t([0-9]+)\t([ATGCN-]+)\t([0-9]+)\t([ATGCN-]+)/),
                dFound ? groomMatch(chunk =~
                        /# Hit table(?:.+\n)+D\t.+\t([0-9]+)\t([ATGCN-]+)\t([0-9]+)\t([ATGCN-]+)/) : null,
                groomMatch(chunk =~
                        /# Hit table(?:.+\n)+J\t.+\t([0-9]+)\t([ATGCN-]+)\t([0-9]+)\t([ATGCN-]+)/)
        ].collect {
            it ? new Alignment(it[0].toInteger() - 1, it[1], it[2].toInteger() - 1, it[3]) : null
        }

        // Deduce CDR/FW markup
        // - CDR1,2 and CDR3(start only) coords
        cdrBounds = [
                groomMatch(chunk =~
                        /# Alignment summary(?:.+\n)+CDR1-IMGT\t([0-9]+)\t([0-9]+)/),
                groomMatch(chunk =~
                        /# Alignment summary(?:.+\n)+CDR2-IMGT\t([0-9]+)\t([0-9]+)/),
                groomMatch(chunk =~
                        /# Alignment summary(?:.+\n)+CDR3-IMGT \(germline\)\t([0-9]+)\t([0-9]+)/)
        ]

        int cdr1Start = -1, cdr1End = -1,
            cdr2Start = -1, cdr2End = -1,
            cdr3Start = -1, cdr3End = -1

        if (cdrBounds[0]) {
            cdr1Start = cdrBounds[0][0].toInteger() - 1
            cdr1End = cdrBounds[0][1].toInteger()
        }

        if (cdrBounds[1]) {
            cdr2Start = cdrBounds[1][0].toInteger() - 1
            cdr2End = cdrBounds[1][1].toInteger()
        }

        // - Find CDR3 end using J reference point manually
        if (cdrBounds[2] && jFound) {
            cdr3Start = cdrBounds[2][0].toInteger() - 4
            cdr3End = jFound ? jRefSearcher.getCdr3End(jSegment, alignments[2]) : -1
        }

        def regionMarkup = new RegionMarkup(cdr1Start, cdr1End, cdr2Start, cdr2End, cdr3Start, cdr3End)
        def hasCdr3 = cdr3Start >= 0, complete = cdr3End >= 0 && hasCdr3

        // - Markup of V/D/J within CDR3
        int vCdr3End = -1, dCdr3Start = -1, dCdr3End = -1, jCdr3Start = -1,
            vDel = -1, dDel5 = -1, dDel3 = -1, jDel = -1

        if (hasCdr3) {
            vCdr3End = alignments[0].qstart + alignments[0].qseq.length() - cdr3Start
            vDel = vSegment.sequence.length() - alignments[0].send

            if (dFound) {
                dCdr3Start = alignments[1].qstart - cdr3Start
                dCdr3End = dCdr3Start + alignments[1].qseq.length()
                dDel5 = alignments[1].sstart
                dDel3 = dSegment.sequence.length() - alignments[1].send
            }

            if (jFound) {
                jCdr3Start = alignments[2].qstart - cdr3Start
                jDel = alignments[2].sstart
            }
        } else {
            noCdr3.incrementAndGet()
        }

        def cdr3Markup = new Cdr3Markup(vCdr3End, dCdr3Start, dCdr3End, jCdr3Start)
        def truncations = new Truncations(vDel, dDel5, dDel3, jDel)

        // Finally, deal with hypermutations
        // offset for converting coordinate in read to coordinate in germline V
        def mutationExtractor = new MutationExtractor(vSegment, alignments[0], regionMarkup)

        if (hasCdr3) {
            if (dFound) {
                mutationExtractor.extractD(dSegment, alignments[1])
            }
            mutationExtractor.extractJ(jSegment, alignments[2])
        }

        return new Mapping(vSegment, dSegment, jSegment,
                regionMarkup, cdr3Markup, truncations,
                rc, complete, hasCdr3, inFrame, noStop, dFound,
                mutationExtractor.mutations)
    }
}