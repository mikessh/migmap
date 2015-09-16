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

import com.antigenomics.higblast.mutation.Mutation
import com.antigenomics.higblast.mutation.MutationType
import com.antigenomics.higblast.mutation.SubRegion

class RegionMarkup {
    final static RegionMarkup DUMMY = new RegionMarkup(-1, -1, -1, -1, -1, -1)

    final int cdr1Start, cdr1End, cdr2Start, cdr2End, cdr3Start, cdr3End // in read, used for extraction

    RegionMarkup(int cdr1Start, int cdr1End,
                 int cdr2Start, int cdr2End,
                 int cdr3Start, int cdr3End) {
        this.cdr1Start = cdr1Start
        this.cdr1End = cdr1End
        this.cdr2Start = cdr2Start
        this.cdr2End = cdr2End
        this.cdr3Start = cdr3Start
        this.cdr3End = cdr3End
    }

    SubRegion getSubRegion(Mutation mutation) {
        getSubRegion(mutation.type == MutationType.Deletion ? mutation.posInRead - 1 : mutation.posInRead)
    }

    SubRegion getSubRegion(int pos) {
        if (pos < cdr1Start) {
            return SubRegion.FR1
        } else if (pos < cdr1End) {
            return SubRegion.CDR1
        } else if (pos < cdr2Start) {
            return SubRegion.FR2
        } else if (pos < cdr2End) {
            return SubRegion.CDR2
        } else if (pos < cdr3Start) {
            return SubRegion.FR3
        } else {
            return SubRegion.CDR3
        }
    }

    boolean isComplete() {
        cdr1Start >= 0 && cdr1End >= 0 &&
                cdr2Start >= 0 && cdr2End >= 0 &&
                cdr3Start >= 0
    }

    static final String OUTPUT_HEADER = "cdr1.start.in.read\tcdr1.end.in.read\t" +
            "cdr2.start.in.read\tcdr2.end.in.read\t" +
            "cdr3.start.in.read\tcdr3.end.in.read"

    @Override
    public String toString() {
        [cdr1Start, cdr1End, cdr2Start, cdr2End, cdr3Start, cdr3End].join("\t")
    }
}
