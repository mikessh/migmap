/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
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
    final static int MIN_CDR3_LEN = 6
    final AtomicInteger total = new AtomicInteger(),
                        noMatch = new AtomicInteger(),
                        vNotFound = new AtomicInteger(),
                        noCdr3 = new AtomicInteger()
    final SegmentDatabase segmentDatabase

    public BlastParser(SegmentDatabase segmentDatabase) {
        this.segmentDatabase = segmentDatabase
    }

    Mapping parse(String chunk) {
        total.incrementAndGet()

        def summary, alignments, cdrBounds

        // Rearrangement summary
        boolean hasD = segmentDatabase.hasD
        if (hasD) {
            summary = groomMatch(chunk =~
                    //         V     D     J  chain stop frame (prod) strand
                    /# V-.+\n(.+)\t(.+)\t(.+)\tV.\t(.+)\t(.+)\t.+\t(.+)\n/)
            if (summary == null) {
                hasD = false
                summary = groomMatch(chunk =~
                        //         V     J    chain   stop frame (prod)  strand
                        /# V-.+\n(.+)\t(.+)\tV.\t(.+)\t(.+)\t.+\t(.+)\n/)
            }
        } else {
            summary = groomMatch(chunk =~
                    //         V     J    chain   stop frame (prod)  strand
                    /# V-.+\n(.+)\t(.+)\tV.\t(.+)\t(.+)\t.+\t(.+)\n/)
        }

        if (summary == null) {
            noMatch.incrementAndGet()
            return null
        }

        def rc = summary[-1] != "+", inFrame = summary[-2] != "Out-of-frame", noStop = summary[-3] != "Yes"

        // Information on segments mapped
        // - Segment names, can be multiple of them
        def vSegmentNames = summary[0],
            jSegmentNames = hasD ? summary[2] : summary[1],
            dSegmentNames = hasD ? summary[1] : BLAST_NA

        def dFound = dSegmentNames != BLAST_NA,
            vFound = vSegmentNames != BLAST_NA,
            jFound = jSegmentNames != BLAST_NA

        if (!vFound) {
            vNotFound.incrementAndGet()
            return null
        }

        def vSegment = segmentDatabase.segments[vSegmentNames.split(",").sort()[0]],
            dSegment = dFound ? segmentDatabase.segments[dSegmentNames.split(",").sort()[0]] : Segment.DUMMY_D,
            jSegment = jFound ? segmentDatabase.segments[jSegmentNames.split(",").sort()[0]] : Segment.DUMMY_J

        // - Alignments for V, D and J segments, remember here and further BLAST coordinates are 1-based
        alignments = [
                groomMatch(chunk =~
                        //                                            qstart     qseq        sstart     sseq
                        /# Hit table(?:.+\n)+V\t$vSegment.regexName\t([0-9]+)\t([ATGCN-]+)\t([0-9]+)\t([ATGCN-]+)/),
                dFound ? groomMatch(chunk =~
                        /# Hit table(?:.+\n)+D\t$dSegment.regexName\t([0-9]+)\t([ATGCN-]+)\t([0-9]+)\t([ATGCN-]+)/) : null,
                jFound ? groomMatch(chunk =~
                        /# Hit table(?:.+\n)+J\t$jSegment.regexName\t([0-9]+)\t([ATGCN-]+)\t([0-9]+)\t([ATGCN-]+)/) : null
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
        if (cdrBounds[2]) {
            cdr3Start = cdrBounds[2][0].toInteger() - 4
        } else if (alignments[0]) {
            // try rescue CDR3
            cdr3Start = RefPointSearcher.getCdr3Start(vSegment, alignments[0])
        }

        if (jFound && cdr3Start >= 0) {
            cdr3End = RefPointSearcher.getCdr3End(jSegment, alignments[2])
        }

        def regionMarkup = new RegionMarkup(cdr1Start, cdr1End, cdr2Start, cdr2End, cdr3Start, cdr3End)
        def hasCdr3 = cdr3Start >= 0 && (cdr3End < 0 || cdr3End - cdr3Start >= MIN_CDR3_LEN),
            complete = cdr3End >= 0 && hasCdr3

        // - Markup of V/D/J within CDR3
        int vCdr3End = -1, dCdr3Start = -1, dCdr3End = -1, jCdr3Start = -1,
            vDel = -1, dDel5 = -1, dDel3 = -1, jDel = -1

        if (hasCdr3) {
            vCdr3End = alignments[0].qend - cdr3Start
            vDel = vSegment.sequence.length() - alignments[0].send

            if (dFound) {
                dCdr3Start = alignments[1].qstart - cdr3Start
                dCdr3End = dCdr3Start + alignments[1].qLength
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
            if (jFound) {
                mutationExtractor.extractJ(jSegment, alignments[2])
            }
        }

        return new Mapping(vSegment, dSegment, jSegment, alignments[0].sstart, alignments[0].qstart,
                regionMarkup, cdr3Markup, truncations,
                rc, complete, hasCdr3, inFrame, noStop, dFound,
                mutationExtractor.mutations)
    }
}