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

import com.antigenomics.higblast.io.Read

class ReadMappingDetailsProvider {
    public static
    final List<String> ALLOWED_FIELDS = Collections.unmodifiableList(["fr1", "cdr1", "fr2", "cdr2", "fr3",
                                                                      "contig"])
    public static ReadMappingDetailsProvider DUMMY = new ReadMappingDetailsProvider([])

    private final String sep
    private final List<String> fields

    ReadMappingDetailsProvider(List<String> fields) {
        this.fields = fields.collect { it.toLowerCase() }
        this.sep = fields.empty ? "" : "\t"
        def badFields = fields.findAll { !ALLOWED_FIELDS.contains(it) }
        if (!badFields.empty) {
            throw new RuntimeException("Bad fields supplied: ${badFields.join(",")}, " +
                    "allowed values: ${ALLOWED_FIELDS.join(",")}.")
        }
    }

    String getHeader() {
        sep + fields.join("\t")
    }

    String getDetails(ReadMapping readMapping) {
        def details = readMapping.mapped ?
                new ReadMappingDetails(readMapping.read, readMapping.mapping) :
                ReadMappingDetails.DUMMY

        sep + fields.collect { details."$it" }.join("\t")
    }

    private static class ReadMappingDetails {
        final String seq
        final RegionMarkup readMarkup, referenceMarkup
        final int vStartInRef

        static final ReadMappingDetails DUMMY = new ReadMappingDetails()

        ReadMappingDetails(Read read, Mapping mapping) {
            this.seq = read.seq
            this.readMarkup = mapping.regionMarkup
            this.vStartInRef = mapping.vStartInRef
            this.referenceMarkup = mapping.vSegment.regionMarkup
        }

        ReadMappingDetails() {
            this.seq = null
            this.vStartInRef = 480011
            this.referenceMarkup = RegionMarkup.DUMMY
            this.readMarkup = RegionMarkup.DUMMY
        }

        String getFr1() {
            vStartInRef < referenceMarkup.cdr1Start ?
                    'N' * vStartInRef + seq.substring(0, readMarkup.cdr1Start) :
                    'N' * referenceMarkup.cdr1Start
        }

        String getCdr1() {
            vStartInRef < referenceMarkup.cdr1End ?
                    'N' * (vStartInRef - referenceMarkup.cdr1Start) + seq.substring(readMarkup.cdr1Start, readMarkup.cdr1End) :
                    'N' * (referenceMarkup.cdr1End - referenceMarkup.cdr1Start)
        }

        String getFr2() {
            vStartInRef < referenceMarkup.cdr2Start ?
                    'N' * (vStartInRef - referenceMarkup.cdr1End) + seq.substring(readMarkup.cdr1End, readMarkup.cdr2Start) :
                    'N' * (referenceMarkup.cdr2Start - referenceMarkup.cdr1End)
        }

        String getCdr2() {
            vStartInRef < referenceMarkup.cdr2End ?
                    'N' * (vStartInRef - referenceMarkup.cdr2Start) + seq.substring(readMarkup.cdr2Start, readMarkup.cdr2End) :
                    'N' * (referenceMarkup.cdr2End - referenceMarkup.cdr2Start)
        }

        String getFr3() {
            vStartInRef < referenceMarkup.cdr3Start ?
                    'N' * (vStartInRef - referenceMarkup.cdr2End) + seq.substring(readMarkup.cdr2End, readMarkup.cdr3Start) :
                    'N' * (referenceMarkup.cdr3Start - referenceMarkup.cdr2End)
        }

        String getContig() {
            fr1 + cdr1 + fr2 + cdr2 + fr3 + downstream
        }

        String getDownstream() {
            readMarkup.cdr3Start >= 0 ? seq.substring(readMarkup.cdr3Start) : ""
        }
    }
}