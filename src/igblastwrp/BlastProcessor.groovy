package igblastwrp

import java.util.regex.Matcher

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
class BlastProcessor {
    final boolean hasD
    final String chain
    final JRefSearcher jRefSearcher

    public BlastProcessor(String chain, JRefSearcher jRefSearcher) {
        this.jRefSearcher = jRefSearcher
        this.hasD = chain =~ /[BH]$/
        this.chain = chain
    }

    private static List<String> groomMatch(Matcher matcher) {
        matcher.size() > 0 ? matcher[0][1..-1] : null//[]
    }

    public Clonotype processChunk(String chunk) {
        def segments

        // Rearrangement summary
        //                                               V     D     J    chain   stop frame prod strand
        //                                               |     |     |      |       |     |    |   |
        segments = hasD ? groomMatch(chunk =~ /# V-.+\n(.+)\t(.+)\t(.+)\tV$chain\t(.+)\t(.+)\t.+\t.+/) :
                //                             V     J    chain   stop frame prod strand
                //                             |     |      |       |     |    |   |
                groomMatch(chunk =~ /# V-.+\n(.+)\t(.+)\tV$chain\t(.+)\t(.+)\t.+\t.+/)


        if (segments == null)
            return // not match for given chain

        def V_SEGM = segments[0], J_SEGM = hasD ? segments[2] : segments[1], D_SEGM = hasD ? segments[1] : "N/A"

        if (V_SEGM == "N/A" || J_SEGM == "N/A")
            return

        def J_SEGM_UNIQ = J_SEGM.split(",")[0]

        def dFound = D_SEGM != "N/A"

        // Hits
        def hits = [
                groomMatch(chunk =~ /# Hit table(?:.+\n)+V\t.+\t([0-9]+)\t([ATGC]+)\t([0-9]+)\t([ATGC]+)\n/),
                dFound ? groomMatch(chunk =~ /# Hit table(?:.+\n)+D\t.+\t([0-9]+)\t([ATGC]+)\t([0-9]+)\t([ATGC]+)\n/) : null,
                groomMatch(chunk =~ /# Hit table(?:.+\n)+J\t.+\t([0-9]+)\t([ATGC]+)\t([0-9]+)\t([ATGC]+)\n/)
        ]

        //println hits

        //# V-(D)-J junction
        // Junction sequence
        // Combine sequence inside
        //def junction = (String) (chunk =~ /# V-.+junction.+\n(.+)\n/)[0][0]
        //def splitJunction = junction.replaceAll(/(?:N\\/A)|[\(\)]+/, "").trim().split("\t")
        //junction = splitJunction[1..<splitJunction.length - 1].join("")

        //println junction

        // CDR coords
        def cdrBounds = [
                groomMatch(chunk =~ /# Alignment summary(?:.+\n)+CDR1-IMGT\t([0-9]+)\t([0-9]+)\t/),
                groomMatch(chunk =~ /# Alignment summary(?:.+\n)+CDR2-IMGT\t([0-9]+)\t([0-9]+)\t/),
                groomMatch(chunk =~ /# Alignment summary(?:.+\n)+CDR3-IMGT \(germline\)\t([0-9]+)\t([0-9]+)\t/)
        ]

        //println cdrBounds

        // Extract CDR3
        // REMEMBER coordinates in BLAST output are 1-based
        int cdr1Start = -1, cdr1End = -1,
            cdr2Start = -1, cdr2End = -1,
            cdr3Start = -1, cdr3End = -1

        if (cdrBounds[0]) {
            cdr1Start = Integer.parseInt(cdrBounds[0][0]) - 1
            cdr1End = Integer.parseInt(cdrBounds[0][1]) - 1
        }

        if (cdrBounds[1]) {
            cdr2Start = Integer.parseInt(cdrBounds[1][0]) - 1
            cdr2End = Integer.parseInt(cdrBounds[1][1]) - 1
        }

        if (cdrBounds[2]) {
            cdr3Start = Integer.parseInt(cdrBounds[2][0]) - 4
            cdr3End = jRefSearcher.getJRefPoint(J_SEGM_UNIQ, Integer.parseInt(hits[2][0]) - 1,
                    hits[2][1], Integer.parseInt(hits[2][2]) - 1, hits[2][3]) + 4
        }

        new Clonotype(V_SEGM, D_SEGM, J_SEGM, cdr1Start, cdr1End, cdr2Start, cdr2End, cdr3Start, cdr3End)
    }
}
