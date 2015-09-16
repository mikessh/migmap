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

package com.antigenomics.higblast.mutation

import com.antigenomics.higblast.genomic.Segment

class Mutation {
    final int startInRead, endInRead
    final String ntFrom, ntTo

    int start, end

    String aaFrom, aaTo // todo

    final MutationType type
    Segment parent
    SubRegion subRegion

    Mutation(MutationType type,
             int start, int end, int startInRead, int endInRead,
             String ntFrom, String ntTo) {
        this.type = type
        this.end = end
        this.startInRead = startInRead
        this.endInRead = endInRead
        this.ntFrom = ntFrom
        this.ntTo = ntTo
        this.start = start
    }

    int getPos() {
        type == MutationType.Insertion ? end : start
    }

    int getPosInRead() {
        type == MutationType.Deletion ? endInRead : startInRead // deletion in ref == insertion in query
    }

    @Override
    public String toString() {
        // insertion occur before "position", so that seq[position, ...] is shifted forward
        type.shortName + pos + ":" +
                (type == MutationType.Insertion ? "" : ntFrom) +
                (type == MutationType.Substitution ? ">" : "") +
                (type == MutationType.Deletion ? "" : ntTo)
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Mutation mutation = (Mutation) o

        start == mutation.start &&
                ntFrom == mutation.ntFrom &&
                ntTo == mutation.ntTo &&
                parent == mutation.parent
    }

    @Override
    int hashCode() {
        int result
        result = ntFrom.hashCode()
        result = 31 * result + ntTo.hashCode()
        result = 31 * result + start
        31 * result + parent.hashCode()
    }
}
