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

package com.antigenomics.migmap.io

import com.antigenomics.migmap.Util
import com.antigenomics.migmap.clonotype.Clonotype
import com.antigenomics.migmap.clonotype.ClonotypeAccumulator
import com.antigenomics.migmap.mapping.ReadMapping
import com.antigenomics.migmap.mapping.ReadMappingDetailsProvider

class ClonotypeOutput implements InputPort<ReadMapping> {
    final ReadMappingDetailsProvider readMappingDetailsProvider
    final ClonotypeAccumulator clonotypeAccumulator
    final PlainTextOutput plainTextOutput

    ClonotypeOutput(PlainTextOutput plainTextOutput = StdOutput.INSTANCE,
                    ReadMappingDetailsProvider readMappingDetailsProvider = ReadMappingDetailsProvider.DUMMY) {
        this.plainTextOutput = plainTextOutput
        this.clonotypeAccumulator = new ClonotypeAccumulator()
        this.readMappingDetailsProvider = readMappingDetailsProvider
        if (plainTextOutput != StdOutput.INSTANCE)
            plainTextOutput.put(Clonotype.OUTPUT_HEADER + readMappingDetailsProvider.header)
    }

    @Override
    void put(ReadMapping readMapping) {
        clonotypeAccumulator.put(readMapping)
    }

    @Override
    void close() {
        Util.report("Sorting ${clonotypeAccumulator.clonotypeMap.size()} clonotype entries, " +
                "writing output.", 2)

        clonotypeAccumulator.clonotypeMap.collect {
            new Clonotype(it.key, it.value, clonotypeAccumulator.total)
        }.sort().each {
            plainTextOutput.put(it.toString() + readMappingDetailsProvider.getDetailsString(it.representativeMapping))
        }
        plainTextOutput.close()

        Util.report("Finished.", 2)
    }
}
