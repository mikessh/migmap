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

import static com.antigenomics.higblast.Util.GAP

class RefPointSearcher {
    final static int MAX_J_FR4_TRUNCATIONS = 5, MAX_V_FR4_TRUNCATIONS = 5

    RefPointSearcher() {
    }

    static int getCdr3End(Segment jSegment, Alignment alignment) {
        convertPositionJ(jSegment.referencePoint + 4,
                alignment.qstart, alignment.qseq, alignment.sstart, alignment.sseq,
                alignment.send)
    }

    static int getCdr3Start(Segment vSegment, Alignment alignment) {
        convertPositionV(vSegment.referencePoint - 3,
                alignment.qstart, alignment.qseq, alignment.sstart, alignment.sseq,
                alignment.send, alignment.qend)
    }

    static int convertPositionJ(int pos, int qstart, String qseq, int sstart, String sseq, int send) {
        if (pos > send) // incomplete, F/W out of scope
            return -1

        def truncations = sstart - pos
        if (truncations >= 0) // conserved residue truncated, alignment starts after F/W
            return truncations >= MAX_J_FR4_TRUNCATIONS ? -1 : qstart

        convertPosition(pos, qstart, qseq, sstart, sseq)
    }

    static int convertPositionV(int pos, int qstart, String qseq, int sstart, String sseq, int send, int qend) {
        if (sstart > pos)  // incomplete, C out of scope
            return -1

        def truncations = pos - send
        if (truncations >= 0) // conserved residue truncated, alignment ends before C
            return truncations >= MAX_V_FR4_TRUNCATIONS ? -1 : qend

        convertPosition(pos, qstart, qseq, sstart, sseq)
    }

    static int convertPosition(int pos, int qstart, String qseq, int sstart, String sseq) {
        int posRel = pos - sstart, posDelta = 0 // normal, count deletions
        for (int i = 0; i < posRel; i++) {
            if (qseq[i] == GAP)
                posDelta--
            else if (sseq[i] == GAP) {
                posDelta--
                posRel++
            }
        }

        posRel + posDelta + qstart
    }
}
