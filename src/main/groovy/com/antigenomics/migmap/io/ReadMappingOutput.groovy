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

import com.antigenomics.migmap.mapping.ReadMapping
import com.antigenomics.migmap.mapping.ReadMappingDetailsProvider
import groovy.transform.CompileStatic

@CompileStatic
class ReadMappingOutput implements InputPort<ReadMapping> {
    final PlainTextOutput plainTextOutput
    final ReadMappingDetailsProvider readMappingDetailsProvider

    ReadMappingOutput(PlainTextOutput plainTextOutput = StdOutput.INSTANCE,
                      ReadMappingDetailsProvider readMappingDetailsProvider = ReadMappingDetailsProvider.DUMMY) {
        this.plainTextOutput = plainTextOutput
        this.readMappingDetailsProvider = readMappingDetailsProvider
        if (plainTextOutput != StdOutput.INSTANCE)
            plainTextOutput.put(ReadMapping.OUTPUT_HEADER + readMappingDetailsProvider.header)
    }

    @Override
    void put(ReadMapping readMapping) {
        if (readMapping.mapped) {
            plainTextOutput.put(readMapping.toString() + readMappingDetailsProvider.getDetailsString(readMapping))
        }
    }

    @Override
    void close() {
        plainTextOutput.close()
    }
}
