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
import com.antigenomics.higblast.genomic.Segment
import com.antigenomics.higblast.genomic.SegmentType
import com.antigenomics.higblast.genomic.VSegment
import com.antigenomics.higblast.mapping.RegionMarkup

import static com.antigenomics.higblast.mutation.MutationType.Deletion
import static com.antigenomics.higblast.mutation.MutationType.None

class MutationExtractor {
    final RegionMarkup regionMarkup

    MutationExtractor(RegionMarkup regionMarkup = null) {
        this.regionMarkup = regionMarkup
    }

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
                if (type != Deletion) {
                    writeMutation()
                    type = Deletion
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

    List<Mutation> extract(Segment segment,
                           Alignment alignment) {
        def mutations = extract(alignment)

        switch (segment.type) {
            case SegmentType.V:
                mutations.each {
                    it.region = segment
                    it.subRegion = deduceSubRegionV(segment as VSegment, it.start)
                }
                break
            case SegmentType.J:
                mutations.each {
                    it.region = segment
                    it.subRegion = (it.start > segment.referencePoint + 3) ? SubRegion.FR4 : SubRegion.CDR3
                    if (regionMarkup) {
                        it.start += regionMarkup.jStart
                        it.end += regionMarkup.jStart
                    }
                }
                break
            case SegmentType.D:
                mutations.each {
                    it.region = segment
                    it.subRegion = SubRegion.CDR3
                    if (regionMarkup) {
                        it.start += regionMarkup.dStart
                        it.end += regionMarkup.dStart
                    }
                }
                break
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
