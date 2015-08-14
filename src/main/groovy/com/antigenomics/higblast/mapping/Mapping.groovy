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

package com.antigenomics.higblast.mapping

import com.antigenomics.higblast.Util
import com.antigenomics.higblast.blast.ClonotypeData
import com.antigenomics.higblast.genomic.DSegment
import com.antigenomics.higblast.genomic.JSegment
import com.antigenomics.higblast.genomic.VSegment
import com.antigenomics.higblast.shm.Hypermutation

class Mapping {
    final List<VSegment> vSegments
    final List<DSegment> dSegments
    final List<JSegment> jSegments
    final boolean rc, complete, hasCdr3, inFrame, noStop
    final Cdr3Markup cdr3Markup
    final RegionMarkup regionMarkup
    String cdr3nt = Util.MY_NA, cdr3aa = Util.MY_NA
    final List<Hypermutation> hypermutations

    Mapping(List<VSegment> vSegments, List<DSegment> dSegments, List<JSegment> jSegments,
            RegionMarkup regionMarkup, Cdr3Markup cdr3Markup,
            boolean rc, boolean complete, boolean hasCdr3, boolean inFrame, boolean noStop,
            List<Hypermutation> hypermutations) {
        // needed to extract CDR3 sequence
        this.regionMarkup = regionMarkup
        this.rc = rc
        
        // main data
        this.vSegments = vSegments
        this.dSegments = dSegments //== Util.BLAST_NA ? Util.MY_NA : dSegment
        this.jSegments = jSegments

        this.cdr3Markup = cdr3Markup
        this.hypermutations = hypermutations

        // misc
        this.complete = complete
        this.hasCdr3 = hasCdr3
        this.inFrame = inFrame
        this.noStop = noStop
    }

    void extractCdr3(String seq) {
        if (hasCdr3) {
            if (rc)
                seq = Util.revCompl(seq)

            this.cdr3nt =
                    complete ? seq.substring(cdr3start, cdr3end) : (seq.substring(cdr3start) + "_")

            this.cdr3aa =
                    complete ? Util.translateCdr(cdr3nt) : (Util.translateLinear(cdr3nt) + "_")
        }
    }

    String getSignature() {
        vSegment + "_" + cdr3nt + "_" + jSegment + "_" + hypermutations.collect { it.signature }.join("_")
    }

    String generateKey(String seq, int level, boolean funcOnly, boolean completeOnly, boolean reportNoCdr3) {
        if ((!reportNoCdr3 && !hasCdr3) || (completeOnly && !complete))
            return null

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

        if (funcOnly && (!inFrame || !noStop))
            return null

        switch (level) {
            case 2:
                return [//vSegment, dSegment, jSegment,
                        cdr1nt, cdr2nt, cdr3nt,
                        cdr1aa, cdr2aa, cdr3aa,
                        inFrame, noStop, complete,
                        hypermutations.join('|')
                ].join('\t')
            case 1:
                boolean inFrame = ![cdr1aa, cdr2aa, cdr3aa].any { it.contains("?") },
                        noStop = ![cdr1aa, cdr2aa, cdr3aa].any { it.contains("*") },
                        complete = ![cdr1aa, cdr2aa, cdr3aa].any { it == Util.MY_NA } &&
                                !cdr3aa.contains("_")

                return [//vSegment, dSegment,jSegment,
                        cdr1nt, cdr2nt, cdr3nt,
                        cdr1aa, cdr2aa, cdr3aa,
                        inFrame, noStop, complete,
                        Util.MY_NA
                ].join('\t')
            default:
                boolean inFrame = !cdr3aa.contains("?"), noStop = !cdr3aa.contains("*"),
                        complete = !cdr3aa.contains("_")

                return [//vSegment, dSegment,jSegment,
                        Util.MY_NA, Util.MY_NA, cdr3nt,
                        Util.MY_NA, Util.MY_NA, cdr3aa,
                        inFrame, noStop, complete,
                        Util.MY_NA
                ].join('\t')
        }
    }

    ClonotypeData appendToData(ClonotypeData clonotypeData, String qual, byte qualThreshold, int nReads, int nEvents, int level) {
        if (rc && qual)
            qual = qual.reverse()

        def filteredHyperm = qual ? hypermutations.findAll {
            ((int) qual.charAt(it.posInRead) - 33) >= qualThreshold
        } : hypermutations

        def cdr1q = level > 0 && cdr1start >= 0 && qual ? qual.substring(cdr1start, cdr1end) : Util.MY_NA,
            cdr2q = level > 0 && cdr2start >= 0 && qual ? qual.substring(cdr2start, cdr2end) : Util.MY_NA,
            cdr3q = cdr3start >= 0 && qual ?
                    (complete ? qual.substring(cdr3start, cdr3end) : qual.substring(cdr3start))
                    : Util.MY_NA

        if (clonotypeData) {
            clonotypeData.append(cdr1q, cdr2q, cdr3q, filteredHyperm, vSegment, dSegment, jSegment, nReads, nEvents)
            return null
        }

        return new ClonotypeData(cdr1q, cdr2q, cdr3q, filteredHyperm, vSegment, dSegment, jSegment, nReads, nEvents, level)
    }

    final static String KEY_HEADER = "cdr1nt\tcdr2nt\tcdr3nt\t" +
            "cdr1aa\tcdr2aa\tcdr3aa\tinFrame\tnoStop\tcomplete"

    final static String HEADER = "v_segment\td_segment\tj_segment\t" +
            "cdr1start\tcdr1end\t" +
            "cdr2start\tcdr2end\t" +
            "cdr3start\tcdr3end\t" +
            "inFrame\tnoStop\tcomplete"

    @Override
    String toString() {
        [cdr3nt, cdr3aa, vSegment, dSegment, jSegment,
         vEnd, dStart, dEnd, jStart,
         hypermutations.collect { it.toString() }].flatten().join("\t")
    }
}
