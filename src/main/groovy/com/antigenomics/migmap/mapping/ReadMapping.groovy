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
import com.antigenomics.migmap.io.Read
import com.antigenomics.migmap.mutation.MutationType

class ReadMapping {
    final Read read
    final Mapping mapping
    final String cdr3nt, cdr3aa
    final byte[] mutationQual, cdrInsertQual
    final boolean canonical, inFrame, noStop

    ReadMapping(Mapping mapping, Read read) {
        this.mapping = mapping
        this.read = read

        if (mapping) {
            if (mapping.rc) {
                // don't forget to reverse complement
                read = read.rc
            }

            def seq = read.seq
            def regionMarkup = mapping.regionMarkup
            
            this.mutationQual = new byte[mapping.mutations.size()]

            mapping.mutations.eachWithIndex { it, i ->
                mutationQual[i] = it.type == MutationType.Substitution ? read.qualAt(it.posInRead) : Util.MAX_QUAL
            }

            def cdrMarkup = mapping.cdr3Markup

            if (mapping.hasCdr3) {
                if (mapping.complete) {
                    this.cdr3nt = seq.substring(regionMarkup.cdr3Start, regionMarkup.cdr3End)
                    this.cdr3aa = Util.translateCdr(cdr3nt)

                    int i = 0

                    if (mapping.dFound) {
                        this.cdrInsertQual = new byte[Math.max(0, cdrMarkup.dStart - cdrMarkup.vEnd) +
                                Math.max(0, cdrMarkup.jStart - cdrMarkup.dEnd)]

                        if (cdrMarkup.vEnd <= cdrMarkup.dStart) {
                            ((regionMarkup.cdr3Start + cdrMarkup.vEnd)..<(regionMarkup.cdr3Start + cdrMarkup.dStart)).each {
                                cdrInsertQual[i++] = read.qualAt(it)
                            }
                        }

                        if (cdrMarkup.dEnd <= cdrMarkup.jStart) {
                            ((regionMarkup.cdr3Start + cdrMarkup.dEnd)..<(regionMarkup.cdr3Start + cdrMarkup.jStart)).each {
                                cdrInsertQual[i++] = read.qualAt(it)
                            }
                        }
                    } else {
                        this.cdrInsertQual = new byte[Math.max(0, cdrMarkup.jStart - cdrMarkup.vEnd)]
                        if (cdrMarkup.vEnd <= cdrMarkup.jStart) {
                            ((regionMarkup.cdr3Start + cdrMarkup.vEnd)..<(regionMarkup.cdr3Start + cdrMarkup.jStart)).each {
                                cdrInsertQual[i++] = read.qualAt(it)
                            }
                        }
                    }
                    this.canonical = Util.isCanonical(cdr3nt)
                } else {
                    this.cdr3nt = seq.substring(regionMarkup.cdr3Start)
                    this.cdr3aa = Util.translateLinear(cdr3nt)
                    this.cdrInsertQual = new byte[cdr3nt.length()]
                    (0..<cdr3nt.length()).each {
                        cdrInsertQual[it] = read.qualAt(regionMarkup.cdr3Start + it)
                    }
                    this.canonical = false
                }

                def inFrame = !cdr3aa.contains("?"), noStop = !cdr3aa.contains("*")

                this.inFrame = mapping.inFrame && inFrame
                this.noStop = mapping.inFrame && noStop
            } else {
                this.cdr3nt = Util.MY_NA
                this.cdr3aa = Util.MY_NA
                this.cdrInsertQual = new byte[0]
                this.canonical = false
                this.inFrame = mapping.inFrame
                this.noStop = mapping.inFrame
            }
        } else {
            this.cdr3nt = null
            this.cdr3aa = null
            this.mutationQual = null
            this.cdrInsertQual = null
            this.canonical = false
            this.inFrame = false
            this.noStop = false
        }
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
    final String OUTPUT_HEADER = "read.header\tcdr3nt\tcdr3aa\tcdr.insert.qual\tmutations.qual\t" + Mapping.OUTPUT_HEADER + "\tcanonical"

    @Override
    String toString() {
        [read.header, cdr3nt, cdr3aa, Util.qualToString(cdrInsertQual), Util.qualToString(mutationQual), mapping.toString(), canonical].join("\t")
    }
}
