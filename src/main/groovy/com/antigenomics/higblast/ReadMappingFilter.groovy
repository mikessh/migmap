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


package com.antigenomics.higblast

import com.antigenomics.higblast.mapping.ReadMapping

import java.util.concurrent.atomic.AtomicLong

class ReadMappingFilter {
    final boolean allowNoCdr3, allowIncomplete, allowNonCoding
    final byte qualityThreshold

    private final AtomicLong totalCounter = new AtomicLong(),
                             mappedCounter = new AtomicLong(),
                             goodCounter = new AtomicLong(),
                             noCdr3Counter = new AtomicLong(),
                             incompleteCounter = new AtomicLong(),
                             nonCodingCounter = new AtomicLong(),
                             passedCounter = new AtomicLong()

    ReadMappingFilter(byte qualityThreshold, boolean allowNoCdr3, boolean allowIncomplete, boolean allowNonCoding) {
        this.qualityThreshold = qualityThreshold
        this.allowNoCdr3 = allowNoCdr3
        this.allowIncomplete = allowIncomplete
        this.allowNonCoding = allowNonCoding
    }

    ReadMappingFilter() {
        this((byte) 25, false, false, true)
    }

    boolean pass(ReadMapping readMapping) {
        totalCounter.incrementAndGet()

        if (readMapping.mapped) {
            mappedCounter.incrementAndGet()

            def mapping = readMapping.mapping

            boolean passQuality = ((readMapping.minCdrInsertQual >= qualityThreshold &&
                    readMapping.minMutationQual >= qualityThreshold) && goodCounter.incrementAndGet() > 0),
                    passCdr3 = (mapping.hasCdr3 || noCdr3Counter.incrementAndGet() < 0 || allowNoCdr3),
                    passIncomplete = (mapping.complete || incompleteCounter.incrementAndGet() < 0 || allowIncomplete),
                    passNonCoding = ((mapping.inFrame && mapping.noStop) || nonCodingCounter.incrementAndGet() < 0 || allowNonCoding)

            return passQuality && passCdr3 && passIncomplete && passNonCoding && passedCounter.incrementAndGet() > 0
        }

        false
    }


    long getTotal() {
        totalCounter.get()
    }

    long getMapped() {
        mappedCounter.get()
    }

    double getMappedRatio() {
        (double) mapped / total
    }

    long getGood() {
        goodCounter.get()
    }

    double getGoodRatio() {
        (double) good / mapped
    }

    long getNoCdr3() {
        noCdr3Counter.get()
    }

    double getNoCdr3Ratio() {
        (double) noCdr3 / mapped
    }

    long getIncomplete() {
        incompleteCounter.get()
    }

    double getIncompleteRatio() {
        (double) incomplete / mapped
    }

    long getNonCoding() {
        nonCodingCounter.get()
    }

    double getNonCodingRatio() {
        (double) nonCoding / mapped
    }

    long getPassed() {
        passedCounter.get()
    }

    double getPassedRatio() {
        (double) passed / total
    }

    static
    final String OUTPUT_HEADER = "total\tpassed\tmapped\tgood\tno.cdr3\tincomplete\tnon.coding"


    @Override
    public String toString() {
        [total, passed, mapped, good, noCdr3, incomplete, nonCoding].join("\t")
    }

    public String toProgressString() {
        "Processed $total reads, of them ${Util.toPercent(mappedRatio)}% mapped " +
                "and ${Util.toPercent(passedRatio)}% passed all filters."
    }
}
