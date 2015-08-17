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

import com.antigenomics.higblast.InputPort
import com.antigenomics.higblast.PlainTextOutput
import com.antigenomics.higblast.mapping.ReadMapping

class ClonotypeOutput implements InputPort<ReadMapping> {
    final ClonotypeAccumulator clonotypeAccumulator = new ClonotypeAccumulator()
    final PlainTextOutput plainTextOutput
    final byte qualThreshold

    ClonotypeOutput(PlainTextOutput plainTextOutput, byte qualThreshold) {
        this.plainTextOutput = plainTextOutput
        this.qualThreshold = qualThreshold
    }

    @Override
    void put(ReadMapping obj) {
        clonotypeAccumulator.put(obj)
    }

    @Override
    void close() {
        clonotypeAccumulator.clonotypeMap.collect {
            new Clonotype(it.key, it.value, clonotypeAccumulator.total)
        }.sort().each {
            if (it.minQual >= qualThreshold) {
                plainTextOutput.put(it.toString())
            }
        }
    }
}
