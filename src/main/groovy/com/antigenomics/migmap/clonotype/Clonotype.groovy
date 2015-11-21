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

package com.antigenomics.migmap.clonotype

import com.antigenomics.migmap.Util
import com.antigenomics.migmap.blast.PSegments
import com.antigenomics.migmap.genomic.Segment
import com.antigenomics.migmap.mapping.Cdr3Markup
import com.antigenomics.migmap.mapping.ReadMapping
import com.antigenomics.migmap.mapping.Truncations
import com.antigenomics.migmap.mutation.Mutation
import com.antigenomics.migmap.mutation.MutationStringifier
import groovy.transform.CompileStatic

@CompileStatic
class Clonotype implements Comparable<Clonotype> {
    final String cdr3nt, cdr3aa
    final Segment vSegment, dSegment, jSegment
    final List<Mutation> mutations
    final long count
    final double freq
    final byte[] cdrInsertQual, mutationQual
    final boolean hasCdr3, complete, inFrame, noStop, canonical
    final Cdr3Markup cdr3Markup
    final Truncations truncations
    final PSegments pSegments
    final ReadMapping representativeMapping

    Clonotype(ClonotypeKey key, ClonotypeData data, long total) {
        this.representativeMapping = key.representativeMapping
        this.cdr3nt = representativeMapping.cdr3nt
        this.cdr3aa = representativeMapping.cdr3aa
        this.vSegment = key.vSegment
        this.dSegment = key.dSegment
        this.jSegment = key.jSegment
        this.mutations = key.mutations
        this.count = data.count.longValue()
        this.freq = (double) count / total

        this.cdrInsertQual = new byte[data.cdrInsertQual.length()]
        def cdrQualCount = data.cdrQualCount.longValue()
        (0..<cdrInsertQual.length).each {
            def qual = (byte) ((double) data.cdrInsertQual.get(it) / cdrQualCount)
            cdrInsertQual[it] = qual
        }
        this.mutationQual = new byte[data.mutationQual.length()]
        (0..<mutationQual.length).each {
            def qual = (byte) ((double) data.mutationQual.get(it) / count)
            mutationQual[it] = qual
        }

        this.cdr3Markup = representativeMapping.mapping.cdr3Markup
        this.truncations = representativeMapping.mapping.truncations
        this.pSegments = representativeMapping.pSegments
        this.hasCdr3 = representativeMapping.mapping.hasCdr3
        this.complete = representativeMapping.mapping.complete

        this.inFrame = representativeMapping.inFrame
        this.noStop = representativeMapping.noStop
        this.canonical = representativeMapping.canonical
    }

    @Override
    int compareTo(Clonotype o) {
        -this.count.compareTo(o.count)
    }

    static
    final String OUTPUT_HEADER = "freq\tcount\tv\td\tj\tcdr3nt\tcdr3aa\t" + MutationStringifier.OUTPUT_HEADER +
            "\tcdr.insert.qual\tmutations.qual\t" + Cdr3Markup.OUTPUT_HEADER + "\t" + Truncations.OUTPUT_HEADER +
            "\t" + PSegments.OUTPUT_HEADER +
            "\thas.cdr3\tin.frame\tno.stop\tcomplete\tcanonical"

    @Override
    String toString() {
        [freq, count,
         vSegment.toString(), dSegment.toString(), jSegment.toString(),
         cdr3nt, cdr3aa,
         MutationStringifier.toString(mutations),
         Util.qualToString(cdrInsertQual), Util.qualToString(mutationQual),
         cdr3Markup, truncations, pSegments,
         hasCdr3, inFrame, noStop, complete, canonical].join("\t")
    }
}
