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

package com.antigenomics.higblast.shm

import com.antigenomics.higblast.Util
import com.antigenomics.higblast.blast.Alignment
import com.antigenomics.higblast.genomic.DSegment
import com.antigenomics.higblast.genomic.JSegment
import com.antigenomics.higblast.genomic.Segment
import com.antigenomics.higblast.genomic.VSegment
import com.antigenomics.higblast.mapping.Cdr3Markup

import static com.antigenomics.higblast.shm.MutationType.None

class MutationExtractor {

    List<Mutation> extract(Alignment alignment) {
        def mutations = new LinkedList<Mutation>()
        int start = -1, end = -1
        int qdelta = 0, sdelta = 0
        def type = None

        def writeMutation = {
            if (type != None) {
                mutations.add(new Mutation(
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
                if (type != MutationType.Substiotution) {
                    writeMutation()
                }
                type = MutationType.Substiotution
                start = i
                end = i + 1
                writeMutation()
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

    List<Mutation> extract(Segment segment,Alignment alignment,
                           Cdr3Markup cdr3Markup // to deduce D frame
                           ) {
        def mutations = new LinkedList<Mutation>()

        int qdelta = 0, sdelta = 0
        for (int i = 0; i < alignment.qseq.length(); i++) {
            char q = alignment.qseq.charAt(i), s = alignment.sseq.charAt(i)

            if (q == Util.GAP) {
                qdelta++
            } else if (s == Util.GAP) {
                sdelta++
            } else if (s != q) {
                int pos = i + sstart - sdelta, posInRead = i + qstart - qdelta,
                    codonPos = pos / 3,
                    codonNtStart = codonPos * 3, codonNtPos = codonNtStart + (pos % 3)

                char aaFrom = '?', aaTo = '?'

                if (codonPos < vSeq.aaSeq.length()) {
                    aaFrom = vSeq.aaSeq.charAt(codonPos)

                    char[] codon = new char[3]
                    for (int j = 0; j < 3; j++) {
                        int jj = codonNtStart + j
                        if (jj == codonNtPos)
                            codon[j] = q
                        else
                            codon[j] = vSeq.ntSeq.charAt(jj)
                    }
                    aaTo = Util.codon2aa(new String(codon))
                }

                def region = null

                if (cdr1start >= 0) {
                    if (posInRead < cdr1start)
                        region = "FW1"
                    else if (posInRead < cdr1end)
                        region = "CDR1"


                    if (cdr2start >= 0) {
                        if (posInRead < cdr2start)
                            region = "FW2"
                        else if (posInRead < cdr2end)
                            region = "CDR2"
                        else
                            region = "FW3"
                    }
                }

                if (!region)
                    region = "NA"

                hypermutations.add(new Hypermutation(pos, posInRead,
                        s, q, aaFrom, aaTo,
                        region))
            }
        }

        hypermutations
    }

    SubRegion deduceSubRegion(Segment segment, int pos) {
        if (segment instanceof VSegment) {
            if (pos < segment.cdr1start) {
                return SubRegion.FR1
            } else if (pos < segment.cdr2end) {
                return SubRegion.CDR1
            } else if (pos < segment.cdr2start) {
                return SubRegion.FR2
            } else if (pos < segment.cdr1end) {
                return SubRegion.CDR2
            } else if (pos < segment.referencePoint - 3) {
                return SubRegion.FR3
            } else {
                return SubRegion.CDR3
            }
        } else if (segment instanceof JSegment) {
            return pos > segment.referencePoint + 3 ? SubRegion.FR4 : SubRegion.CDR3
        } else if (segment instanceof DSegment) {
            return SubRegion.CDR3
        }

        null
    }
}
