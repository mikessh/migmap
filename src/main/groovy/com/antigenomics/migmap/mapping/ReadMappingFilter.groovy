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
import com.antigenomics.migmap.io.DummyInputPort
import com.antigenomics.migmap.io.InputPort
import com.antigenomics.migmap.io.Read

import java.util.concurrent.atomic.AtomicLong

class ReadMappingFilter {
    final boolean allowNoCdr3, allowIncomplete, allowNonCoding, allowNonCanonical
    final byte qualityThreshold
    final InputPort<Read> unmappedInputPort

    final static ReadMappingFilter createDummy() { new ReadMappingFilter((byte) 0, true, true, true, true) }

    private final AtomicLong totalCounter = new AtomicLong(),
                             mappedCounter = new AtomicLong(),
                             goodCounter = new AtomicLong(),
                             noCdr3Counter = new AtomicLong(),
                             incompleteCounter = new AtomicLong(),
                             nonCodingCounter = new AtomicLong(),
                             nonCanonicalCounter = new AtomicLong(),
                             passedCounter = new AtomicLong()

    ReadMappingFilter(byte qualityThreshold, boolean allowNoCdr3, boolean allowIncomplete, boolean allowNonCoding,
                      boolean allowNonCanonical,
                      InputPort<Read> unmappedInputPort = DummyInputPort.INSTANCE) {
        this.qualityThreshold = qualityThreshold
        this.allowNoCdr3 = allowNoCdr3
        this.allowIncomplete = allowIncomplete
        this.allowNonCoding = allowNonCoding
        this.allowNonCanonical = allowNonCanonical
        this.unmappedInputPort = unmappedInputPort
    }

    ReadMappingFilter() {
        this((byte) 25, false, false, true, false)
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
                    passNonCoding = ((mapping.inFrame && mapping.noStop) || nonCodingCounter.incrementAndGet() < 0 || allowNonCoding),
                    passNonCanonical = (readMapping.canonical || nonCanonicalCounter.incrementAndGet() < 0 || allowNonCanonical)

            return passQuality && passCdr3 && passIncomplete && passNonCoding && passNonCanonical && passedCounter.incrementAndGet() > 0
        } else {
            unmappedInputPort.put(readMapping.read)
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

    long getNonCanonical() {
        nonCanonicalCounter.get()
    }

    double getNonCanonicalRatio() {
        (double) nonCanonical / mapped
    }

    long getPassed() {
        passedCounter.get()
    }

    double getPassedRatio() {
        (double) passed / mapped
    }

    static
    final String OUTPUT_HEADER = "total\tpassed\tmapped\tgood.of.mapped\t" +
            "no.cdr3.of.mapped\tincomplete.of.mapped\tnon.coding.of.mapped\tnon.canonical.of.mapped"


    @Override
    public String toString() {
        [total, passed, mapped, good, noCdr3, incomplete, nonCoding, nonCanonical].join("\t")
    }

    public String toProgressString() {
        "Processed $total reads, of them ${Util.toPercent(mappedRatio)}% mapped. " +
                "Among mapped reads ${Util.toPercent(goodRatio)}% passed quality filter " +
                "and ${Util.toPercent(passedRatio)}% passed all filters."
    }
}
