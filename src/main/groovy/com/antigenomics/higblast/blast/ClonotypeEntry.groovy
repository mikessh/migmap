package com.antigenomics.higblast.blast

import com.antigenomics.higblast.shm.Hypermutation

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
class ClonotypeEntry {
    final String vSegment, dSegment, jSegment,
                 cdr1nt, cdr2nt, cdr3nt,
                 cdr1aa, cdr2aa, cdr3aa
    final List<Hypermutation> hypermutations
    final boolean inFrame, noStop, complete

    ClonotypeEntry(String vSegment, String dSegment, String jSegment,
                   String cdr1nt, String cdr2nt, String cdr3nt,
                   String cdr1aa, String cdr2aa, String cdr3aa,
                   boolean inFrame, boolean noStop, boolean complete,
                   List<Hypermutation> hypermutations) {
        this.vSegment = vSegment
        this.dSegment = dSegment
        this.jSegment = jSegment
        this.cdr1nt = cdr1nt
        this.cdr2nt = cdr2nt
        this.cdr3nt = cdr3nt
        this.cdr1aa = cdr1aa
        this.cdr2aa = cdr2aa
        this.cdr3aa = cdr3aa
        this.inFrame = inFrame
        this.noStop = noStop
        this.complete = complete
        this.hypermutations = hypermutations
    }

    final static String HEADER = "v_segment\td_segment\tj_segment\t" +
            "cdr1nt\tcdr2nt\tcdr3nt\t" +
            "cdr1aa\tcdr2aa\tcdr3aa\t" +
            "inFrame\tnoStop\tcomplete"

    @Override
    String toString() {
        [vSegment, dSegment, jSegment,
         cdr1nt, cdr2nt, cdr3nt,
         cdr1aa, cdr2aa, cdr3aa,
         inFrame, noStop, complete].join("\t")
    }

    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        ClonotypeEntry that = (ClonotypeEntry) o

        if (cdr1nt != that.cdr1nt) return false
        if (cdr2nt != that.cdr2nt) return false
        if (cdr3nt != that.cdr3nt) return false
        if (dSegment != that.dSegment) return false
        if (jSegment != that.jSegment) return false
        if (vSegment != that.vSegment) return false

        true
    }

    int hashCode() {
        int result
        result = vSegment.hashCode()
        result = 31 * result + dSegment.hashCode()
        result = 31 * result + jSegment.hashCode()
        result = 31 * result + cdr1nt.hashCode()
        result = 31 * result + cdr2nt.hashCode()
        31 * result + cdr3nt.hashCode()
    }
}
