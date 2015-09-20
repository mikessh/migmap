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

package com.antigenomics.migmap.clonotype

import com.antigenomics.migmap.io.InputPort
import com.antigenomics.migmap.mapping.ReadMapping

import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.atomic.AtomicLong

class ClonotypeAccumulator implements InputPort<ReadMapping> {
    final AtomicLong totalCounter = new AtomicLong()
    final ConcurrentHashMap<ClonotypeKey, ClonotypeData> clonotypeMap = new ConcurrentHashMap<>()

    ClonotypeAccumulator() {
    }

    @Override
    void put(ReadMapping readMapping) {
        totalCounter.incrementAndGet()
        ClonotypeData clonotypeData = clonotypeMap.putIfAbsent(new ClonotypeKey(readMapping),
                new ClonotypeData(readMapping))
        if (clonotypeData) {
            clonotypeData.update(readMapping)
        }
    }

    @Override
    void close() {

    }

    long getTotal() {
        totalCounter.get()
    }
}
