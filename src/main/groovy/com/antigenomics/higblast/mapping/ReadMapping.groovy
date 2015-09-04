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
import com.antigenomics.higblast.io.Read
import com.antigenomics.higblast.mutation.MutationType

class ReadMapping {
    final Read read
    final Mapping mapping
    final String cdr3nt
    final byte[] mutationQual, cdrInsertQual
    final boolean canonical

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
            this.cdr3nt = mapping.hasCdr3 ?
                    (mapping.complete ?
                            seq.substring(regionMarkup.cdr3Start, regionMarkup.cdr3End) :
                            seq.substring(regionMarkup.cdr3Start))
                    : Util.MY_NA

            this.mutationQual = new byte[mapping.mutations.size()]

            mapping.mutations.eachWithIndex { it, i ->
                mutationQual[i] = it.type == MutationType.Substitution ? read.qualAt(it.posInRead) : Util.MAX_QUAL
            }

            def cdrMarkup = mapping.cdr3Markup

            if (mapping.hasCdr3) {
                if (mapping.complete) {
                    int i = 0

                    if (mapping.hasD) {
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
                    this.cdrInsertQual = new byte[cdr3nt.length()]
                    (0..<cdr3nt.length()).each {
                        cdrInsertQual[it] = read.qualAt(regionMarkup.cdr3Start + it)
                    }
                    this.canonical = false
                }
            } else {
                this.cdrInsertQual = new byte[0]
                this.canonical = false
            }
        } else {
            this.cdr3nt = null
            this.mutationQual = null
            this.cdrInsertQual = null
            this.canonical = false
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
    final String OUTPUT_HEADER = "read.header\tcdr3nt\tcdr.insert.qual\tmutations.qual\t" + Mapping.OUTPUT_HEADER + "\tcanonical"

    @Override
    String toString() {
        [read.header, cdr3nt, Util.qualToString(cdrInsertQual), Util.qualToString(mutationQual), mapping.toString(), canonical].join("\t")
    }
}
