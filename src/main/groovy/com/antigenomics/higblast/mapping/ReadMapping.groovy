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

class ReadMapping {
    final Mapping mapping
    final String header, cdr3nt
    final byte[] mutationQual, cdrInsertQual

    ReadMapping(Mapping mapping, Read read) {
        this.mapping = mapping
        this.header = read.header

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
                mutationQual[i] = read.qualAt(it.posInRead)
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
                                cdrInsertQual[i++] = read.qualAt(it) // todo: exception here
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
                } else {
                    this.cdrInsertQual = new byte[cdr3nt.length()]
                    (0..<cdr3nt.length()).each {
                        cdrInsertQual[it] = read.qualAt(regionMarkup.cdr3Start + it)
                    }
                }
            } else {
                this.cdrInsertQual = new byte[0]
            }

        } else {
            this.cdr3nt = null
            this.mutationQual = null
            this.cdrInsertQual = null
        }
    }

    @Override
    String toString() {
        header
    }
}
