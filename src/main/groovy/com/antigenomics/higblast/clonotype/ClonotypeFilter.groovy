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

package com.antigenomics.higblast.clonotype

import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicLong

class ClonotypeFilter {
    final boolean allowNoCdr3, allowIncomplete, allowNonCoding
    private final AtomicLong totalCounter = new AtomicLong(),
                                noCdr3Counter = new AtomicLong(),
                                incompleteCounter = new AtomicLong(),
                                nonCodingCounter = new AtomicLong(),
                                passedCounter = new AtomicLong()

    ClonotypeFilter(boolean allowNoCdr3, boolean allowIncomplete, boolean allowNonCoding) {
        this.allowNoCdr3 = allowNoCdr3
        this.allowIncomplete = allowIncomplete
        this.allowNonCoding = allowNonCoding
    }

    ClonotypeFilter() {
        this(false, false, true)
    }

    boolean pass(Clonotype clonotype) {
        totalCounter.incrementAndGet()
        boolean passCdr3 = (clonotype.hasCdr3 || noCdr3Counter.incrementAndGet() < 0 || allowNoCdr3),
                passIncomplete = (clonotype.complete || incompleteCounter.incrementAndGet() < 0 || allowIncomplete),
                passNonCoding = ((clonotype.inFrame && clonotype.noStop) || nonCodingCounter.incrementAndGet() < 0 || allowNonCoding)
        passCdr3 && passIncomplete && passNonCoding && passedCounter.incrementAndGet() > 0
    }

    long getTotal() {
        totalCounter.get()
    }

    long getNoCdr3() {
        noCdr3Counter.get()
    }

    long getIncomplete() {
        incompleteCounter.get()
    }

    long getNonCoding() {
        nonCodingCounter.get()
    }

    long getPassed() {
        passedCounter.get()
    }

    static
    final String OUTPUT_HEADER = "total\tpassed\tno.cdr3\tincomplete\tnon.coding"


    @Override
    public String toString() {
        [total, passed, noCdr3, incomplete, nonCoding].join("\t")
    }
}
