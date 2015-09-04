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

import static com.antigenomics.higblast.Util.GAP

class JRefSearcher {

    JRefSearcher() {
    }

    static int getCdr3End(Segment jSegment, Alignment alignment) {
        getCdr3End(jSegment.referencePoint, alignment.qstart, alignment.qseq, alignment.sstart, alignment.sseq)
    }

    static int getCdr3End(int jRef, int qstart, String qseq, int sstart, String sseq) {
        if (jRef + 4 > sseq.replaceAll("-", "").length() + sstart) // incomplete, F/W out of scope
            return -1

        if (sstart >= jRef + 4) // conserved residue truncated, alignment starts after F/W
            return qstart

        int jRefRel = jRef - sstart, jRefDelta = 0 // normal, count deletions
        for (int i = 0; i < jRefRel; i++) {
            if (qseq[i] == GAP)
                jRefDelta--
            else if (sseq[i] == GAP) {
                jRefDelta--
                jRefRel++
            }
        }

        jRefRel + jRefDelta + qstart + 4
    }
}
