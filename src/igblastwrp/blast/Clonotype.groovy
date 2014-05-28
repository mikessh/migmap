package igblastwrp.blast

import igblastwrp.Util
import igblastwrp.shm.Hypermutation

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
class Clonotype {
    final String vSegment, dSegment, jSegment
    final int cdr1start, cdr1end, cdr2start, cdr2end, cdr3start, cdr3end
    final boolean rc, inFrame, noStop, complete
    final List<Hypermutation> hypermutations

    Clonotype(String vSegment, String dSegment, String jSegment,
              int cdr1start, int cdr1end, int cdr2start, int cdr2end, int cdr3start, int cdr3end,
              boolean rc, boolean inFrame, boolean noStop, boolean complete,
              List<Hypermutation> hypermutations) {
        this.vSegment = vSegment
        this.dSegment = dSegment == Util.BLAST_NA ? Util.MY_NA : dSegment
        this.jSegment = jSegment
        this.cdr1start = cdr1start
        this.cdr1end = cdr1end
        this.cdr2start = cdr2start
        this.cdr2end = cdr2end
        this.cdr3start = cdr3start
        this.cdr3end = cdr3end
        this.rc = rc
        this.inFrame = inFrame
        this.noStop = noStop
        this.complete = complete
        this.hypermutations = hypermutations
    }

    ClonotypeEntry generateEntry(String seq) {
        if (rc)
            seq = Util.revCompl(seq)

        def cdr1nt = cdr1start >= 0 ? seq.substring(cdr1start, cdr1end) : Util.MY_NA,
            cdr2nt = cdr2start >= 0 ? seq.substring(cdr2start, cdr2end) : Util.MY_NA,
            cdr3nt = cdr3start >= 0 ?
                    (complete ? seq.substring(cdr3start, cdr3end) : seq.substring(cdr3start) + "_")
                    : Util.MY_NA

        def cdr1aa = cdr1start >= 0 ? Util.translateLinear(cdr1nt) : Util.MY_NA,
            cdr2aa = cdr2start >= 0 ? Util.translateLinear(cdr2nt) : Util.MY_NA,
            cdr3aa = cdr3start >= 0 ?
                    (complete ? Util.translateCdr(cdr3nt) : (Util.translateLinear(cdr3nt) + "_"))
                    : Util.MY_NA

        new ClonotypeEntry(vSegment, dSegment, jSegment,
                cdr1nt, cdr2nt, cdr3nt,
                cdr1aa, cdr2aa, cdr3aa,
                inFrame, noStop, complete
        )
    }

    ClonotypeData appendToData(ClonotypeData clonotypeData, String qual, byte qualThreshold) {
        if (rc && qual)
            qual = qual.reverse()

        def filteredHyperm = qual ? hypermutations.findAll {
            ((int)qual.charAt(it.posInRead) - 33) >= qualThreshold
        } : hypermutations

        def cdr1q = cdr1start >= 0 && qual ? qual.substring(cdr1start, cdr1end) : Util.MY_NA,
            cdr2q = cdr2start >= 0 && qual ? qual.substring(cdr2start, cdr2end) : Util.MY_NA,
            cdr3q = cdr3start >= 0 && qual ?
                    (complete ? qual.substring(cdr3start, cdr3end) : qual.substring(cdr3start))
                    : Util.MY_NA

        if (clonotypeData) {
            clonotypeData.append(cdr1q, cdr2q, cdr3q, filteredHyperm)
            return null
        }

        return new ClonotypeData(cdr1q, cdr2q, cdr3q, filteredHyperm)
    }

    boolean isFunctional() {
        inFrame && noStop
    }

    final static String HEADER = "v_segment\td_segment\tj_segment\t" +
            "cdr1start\tcdr1end\t" +
            "cdr2start\tcdr2end\t" +
            "cdr3start\tcdr3end\t" +
            "inFrame\tnoStop\tcomplete"


    @Override
    String toString() {
        [vSegment, dSegment, jSegment,
         cdr1start, cdr1end,
         cdr2start, cdr2end,
         cdr3start, cdr3end,
         inFrame, noStop, complete].join("\t")
    }
}
