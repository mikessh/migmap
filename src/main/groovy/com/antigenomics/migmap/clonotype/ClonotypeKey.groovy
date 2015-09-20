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

import com.antigenomics.migmap.genomic.Segment
import com.antigenomics.migmap.mapping.ReadMapping
import com.antigenomics.migmap.mutation.Mutation

class ClonotypeKey {
    final String cdr3nt
    protected final ReadMapping representativeMapping

    ClonotypeKey(ReadMapping readMapping) {
        this.cdr3nt = readMapping.cdr3nt
        this.representativeMapping = readMapping
    }

    Segment getvSegment() {
        representativeMapping.mapping.vSegment
    }

    Segment getdSegment() {
        representativeMapping.mapping.dSegment
    }

    Segment getjSegment() {
        representativeMapping.mapping.jSegment
    }

    List<Mutation> getMutations() {
        representativeMapping.mapping.mutations
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        ClonotypeKey that = (ClonotypeKey) o

        cdr3nt == that.cdr3nt &&
                vSegment == that.vSegment &&
                jSegment == that.jSegment &&
                mutations == that.mutations
    }

    @Override
    int hashCode() {
        int result
        result = cdr3nt.hashCode()
        result = 31 * result + vSegment.hashCode()
        result = 31 * result + jSegment.hashCode()
        31 * result + mutations.hashCode()
    }
}
