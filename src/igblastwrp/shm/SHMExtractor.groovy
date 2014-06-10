package igblastwrp.shm

import igblastwrp.Util
import igblastwrp.io.FastaReader

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
class SHMExtractor {
    def vSegmentSeqMap = new HashMap<String, VSegmentData>()

    SHMExtractor(String vSegmentFasta) {
        def reader = new FastaReader(vSegmentFasta)

        def read
        while ((read = reader.next()) != null) {
            def vSegment = read.header.substring(1).trim() //trim >
            vSegmentSeqMap.put(vSegment, new VSegmentData(Util.translateLinear(read.seq), read.seq))
        }
    }

    List<Hypermutation> extract(String vSegment, int qstart, String qseq, int sstart, String sseq,
                                int cdr1start, int cdr1end, int cdr2start, int cdr2end) {
        def hypermutations = new LinkedList<Hypermutation>()

        def vSeq = vSegmentSeqMap[vSegment]

        int qdelta = 0, sdelta = 0
        for (int i = 0; i < qseq.length(); i++) {
            char q = qseq.charAt(i), s = sseq.charAt(i)

            if (q == '-') {
                qdelta++
            } else if (s == '-') {
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

                def region = "FW1"
                if (cdr1start >= 0 && posInRead >= cdr1start && posInRead < cdr1end)
                    region = "CDR1"
                else if (cdr2start >= 0) {
                    if (posInRead < cdr2start)
                        region = "FW2"
                    else if (posInRead < cdr2end)
                        region = "CDR2"
                    else
                        region = "FW3"
                }

                hypermutations.add(new Hypermutation(pos, posInRead,
                        s, q, aaFrom, aaTo,
                        region))
            }
        }

        hypermutations
    }
}
