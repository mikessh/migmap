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

package com.antigenomics.migmap.mapping

import com.antigenomics.migmap.Util
import com.antigenomics.migmap.blast.PSegments
import com.antigenomics.migmap.io.Read
import groovy.transform.CompileStatic

@CompileStatic
class ReadMapping {
    final Read read
    final Mapping mapping
    final String cdr3nt, cdr3aa
    final byte[] mutationQual, cdrInsertQual
    final boolean canonical, inFrame, noStop
    final PSegments pSegments

    ReadMapping(Read read, Mapping mapping, String cdr3nt, String cdr3aa,
                byte[] mutationQual, byte[] cdrInsertQual, 
                boolean canonical, boolean inFrame, boolean noStop, 
                PSegments pSegments) {
        this.read = read
        this.mapping = mapping
        this.cdr3nt = cdr3nt
        this.cdr3aa = cdr3aa
        this.mutationQual = mutationQual
        this.cdrInsertQual = cdrInsertQual
        this.canonical = canonical
        this.inFrame = inFrame
        this.noStop = noStop
        this.pSegments = pSegments
    }

    byte getMinCdrInsertQual() {
        Util.minQual(cdrInsertQual)
    }

    byte getMinMutationQual() {
        Util.minQual(mutationQual)
    }

    boolean isMapped() {
        mapping != null
    }

    static
    final String OUTPUT_HEADER = "read.header\tcdr3nt\tcdr3aa\tcdr.insert.qual\tmutations.qual\t" +
            "$Mapping.OUTPUT_HEADER\t$PSegments.OUTPUT_HEADER\tcanonical"

    @Override
    String toString() {
        [read.header, cdr3nt, cdr3aa, Util.qualToString(cdrInsertQual), Util.qualToString(mutationQual),
         mapping.toString(), pSegments.toString(), canonical].join("\t")
    }
}
