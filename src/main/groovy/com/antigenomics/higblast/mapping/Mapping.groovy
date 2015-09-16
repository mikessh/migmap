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

package com.antigenomics.higblast.mapping

import com.antigenomics.higblast.genomic.Segment
import com.antigenomics.higblast.mutation.Mutation
import com.antigenomics.higblast.mutation.MutationStringifier

class Mapping {
    final Segment vSegment, dSegment, jSegment
    final RegionMarkup regionMarkup
    final Cdr3Markup cdr3Markup
    final Truncations truncations
    final List<Mutation> mutations
    final int vStartInRef, vStartInQuery
    final boolean rc, complete, hasCdr3, inFrame, noStop, hasD

    Mapping(Segment vSegment, Segment dSegment, Segment jSegment, int vStartInRef, int vStartInQuery,
            RegionMarkup regionMarkup, Cdr3Markup cdr3Markup, Truncations truncations,
            boolean rc, boolean complete, boolean hasCdr3, boolean inFrame, boolean noStop, boolean hasD,
            List<Mutation> mutations) {
        // needed to extract CDR3 sequence
        this.vStartInRef = vStartInRef
        this.vStartInQuery = vStartInQuery
        this.regionMarkup = regionMarkup
        this.rc = rc

        // main data
        this.vSegment = vSegment
        this.dSegment = dSegment
        this.jSegment = jSegment

        this.cdr3Markup = cdr3Markup
        this.truncations = truncations

        this.mutations = mutations

        // misc
        this.complete = complete
        this.hasCdr3 = hasCdr3
        this.inFrame = inFrame
        this.noStop = noStop
        this.hasD = hasD
    }

    static final String OUTPUT_HEADER = "v.segment\td.segment\tj.segment\t" +
            RegionMarkup.OUTPUT_HEADER + "\t" +
            Cdr3Markup.OUTPUT_HEADER + "\t" +
            Truncations.OUTPUT_HEADER + "\t" +
            MutationStringifier.OUTPUT_HEADER + "\t" +
            "rc\tcomplete\thas.cdr3\tin.frame\tno.stop"


    @Override
    public String toString() {
        [vSegment.toString(), dSegment.toString(), jSegment.toString(),
         regionMarkup.toString(),
         cdr3Markup.toString(),
         truncations.toString(),
         MutationStringifier.toString(mutations),
         rc, complete, hasCdr3, inFrame, noStop].join("\t")
    }
}
