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

package com.antigenomics.higblast.mutation

import com.antigenomics.higblast.Util
import com.antigenomics.higblast.blast.Alignment
import com.antigenomics.higblast.genomic.JSegment
import com.antigenomics.higblast.genomic.Segment
import com.antigenomics.higblast.genomic.VSegment

import static com.antigenomics.higblast.mutation.MutationType.None

class MutationExtractor {

    static List<Mutation> extract(Alignment alignment) {
        def mutations = new LinkedList<Mutation>()
        int start = -1, end = -1
        int qdelta = 0, sdelta = 0
        def type = None

        def writeMutation = {
            if (type != None) {
                mutations.add(new Mutation(
                        type,
                        start + alignment.sstart - sdelta, end + alignment.sstart - sdelta,
                        start + alignment.qstart - qdelta, end + alignment.qstart - qdelta,
                        alignment.sseq[start..<end], alignment.qseq[start..<end])
                )
            }
        }

        for (int i = 0; i < alignment.qseq.length(); i++) {
            char q = alignment.qseq.charAt(i), s = alignment.sseq.charAt(i)

            if (q == Util.GAP) {
                if (type != MutationType.Deletion) {
                    writeMutation()
                    type = MutationType.Deletion
                    start = i
                    end = i + 1
                } else {
                    end++
                }

                qdelta++
            } else if (s == Util.GAP) {
                if (type != MutationType.Insertion) {
                    writeMutation()
                    type = MutationType.Insertion
                    start = i
                    end = i + 1
                } else {
                    end++
                }

                sdelta++
            } else if (s != q) {
                // We store substitutions one-by-one
                // This is to facilitate post analysis (perhaps should be changed)
                if (type != MutationType.Substitution) {
                    writeMutation()
                }
                type = MutationType.Substitution
                start = i
                end = i + 1
                writeMutation()
                type = None
            } else {
                // no mutation here, flush
                writeMutation()
                type = None
            }
        }

        // flush last
        writeMutation()

        mutations
    }

    static List<Mutation> extract(Segment segment,
                                  Alignment alignment) {

        // todo: cdr3start >= 0 !!!!
        def mutations = extract(alignment)

        if (segment instanceof VSegment) {
            mutations.each {
                it.region = segment
                it.subRegion = deduceSubRegionV(segment as VSegment, it.start)
            }
        } else if (segment instanceof JSegment) {
            mutations.each {
                it.region = segment
                it.subRegion = (it.start > segment.referencePoint + 3) ? SubRegion.FR4 : SubRegion.CDR3
            }
        } else {
            mutations.each {
                it.region = segment
                it.subRegion = SubRegion.CDR3
            }
        }

        /*
        int frame = regionMarkup.cdr3Start % 3
        mutations.each {
            int frameShift = (it.startInRead - frame) % 3
            int start = it.startInRead - frameShift,
                end = start + (it.endInRead - it.startInRead) / 3
        }*/

        mutations
    }

    static SubRegion deduceSubRegionV(VSegment segment, int pos) {
        if (pos < segment.cdr1start) {
            return SubRegion.FR1
        } else if (pos < segment.cdr1end) {
            return SubRegion.CDR1
        } else if (pos < segment.cdr2start) {
            return SubRegion.FR2
        } else if (pos < segment.cdr2end) {
            return SubRegion.CDR2
        } else if (pos < segment.referencePoint - 3) {
            return SubRegion.FR3
        } else {
            return SubRegion.CDR3
        }
    }
}
