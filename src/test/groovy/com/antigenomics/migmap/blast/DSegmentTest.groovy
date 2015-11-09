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

package com.antigenomics.migmap.blast

import com.antigenomics.migmap.genomic.Segment
import com.antigenomics.migmap.genomic.SegmentDatabase
import org.junit.AfterClass
import org.junit.Ignore
import org.junit.Test

import static com.antigenomics.migmap.blast.BlastTestUtil.*

class DSegmentTest {
    final BlastInstance blastInstance

    DSegmentTest() {
        blastInstance = new BlastInstanceFactory("data/", "human", ["IGH"], true, false).create()
    }

    @Test
    void igblastOnlyTest() {
        def seq = "GACACGGCCGTATATTTCTGTGCGAGGGATCAGGGGCCACGAGACCACCCCAGTTCATCGTTTTGGGGCCAGGGAACCCTGGTCACCGTC"
        def read = toRead(seq)
        def chunk = toChunk(read)
        def mapping = parser.parse(chunk)
        def readMapping = blastInstance.createReadMapping(mapping, read)

        assert readMapping.mapping.dSegment.name == "IGHD3-3*01"
        assert readMapping.mapping.cdr3Markup.dStart == 60 - 18
        assert readMapping.mapping.cdr3Markup.dEnd == 66 - 18
    }

    @Ignore
    @Test
    void dRecoveryTest() {
        //         ................A........--------------------------------------------
        //         ---------------------------------------------------------------------
        //         ---------------------------------------------........................
        //         ---------------------------------------------------------------------
        //         ---------------------------------------------------------------------
        //         000000000011111111112222222222333333333344444444445555555555666666666
        //         012345678901234567890123456789012345678901234567890123456789012345678
        def seq = "GACACGGCCGTATATTTCTGTGCGAGGGATCAGGGGCCACGAGCCGGCCAGGGAACCCTGGTCACCGTC"
        def read = toRead(seq)
        def chunk = toChunk(read)
        def mapping = parser.parse(chunk)
        def readMapping = blastInstance.createReadMapping(mapping, read)

        assert mapping.dSegment == Segment.DUMMY_D
        assert readMapping.mapping.dSegment.name == "IGHD2-15*01"
        assert readMapping.mapping.cdr3Markup.dStart == 8
        assert readMapping.mapping.cdr3Markup.dEnd == 12
    }

    @AfterClass
    static void tearDown() {
        SegmentDatabase.clearTemporaryFiles()
    }
}
