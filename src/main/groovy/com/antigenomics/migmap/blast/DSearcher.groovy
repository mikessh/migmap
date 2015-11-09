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

package com.antigenomics.migmap.blast

import com.antigenomics.migmap.genomic.Segment
import groovy.transform.CompileStatic

@CompileStatic
class DSearcher {
    final int minHitSize, maxHitSize
    final List<DHitTracker> hitTrackers

    DSearcher(List<Segment> segments) {
        this(segments, 2, 4) // IgBlast has a minimum exact match size of 5
    }

    DSearcher(List<Segment> segments, int minHitSize, int maxHitSize) {
        this.hitTrackers = new LinkedList<>()
        for (Segment segment : segments) {
            hitTrackers.add(new DHitTracker(segment, minHitSize, maxHitSize))
        }
        this.minHitSize = minHitSize
        this.maxHitSize = maxHitSize
    }

    AuxiliaryDMapping scan(String cdr3seq, int vEndRel, int jStartRel) {
        int regionSize = jStartRel - vEndRel - 1
        if (regionSize < minHitSize)
            return null

        def subSeq = cdr3seq.substring(vEndRel + 1, jStartRel)
        def hitSegments = new ArrayList<>()

        // iterate from largest window to smallest one
        for (int i = maxHitSize; i >= minHitSize; i--) {
            // sliding window scan
            for (int j = 0; j < subSeq.length() - i; j++) {
                String kmer = subSeq.substring(j, j + i)

                boolean hasHit = false
                for (DHitTracker hitTracker : hitTrackers) {
                    if (hitTracker.hasHit(kmer)) {
                        hasHit = true
                        hitSegments.add(hitTracker.getSegment())
                    }
                }

                if (hasHit) {
                    int from = vEndRel + 1 + j;
                    return new AuxiliaryDMapping(from,
                            from + i,
                            // binomial P-value, uniform base distribution
                            1.0 - Math.pow(1.0 - Math.pow(0.25, i), regionSize - i + 1),
                            hitSegments
                    )
                }
            }
        }

        null
    }
}
