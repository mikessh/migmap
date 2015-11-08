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

import com.antigenomics.migmap.mapping.Cdr3Markup
import org.junit.Test

class PSegmentTest {
    @Test
    void noPTest() {
        //                    V1     D0   D1           J0
        //                    |      |    |            |
        //             000000000011111111112222222222333333333344444
        //             012345678901234567890123456789012345678901234
        def cdr3seq = "TGTGCGAGGTGGCTTGGGGAAGACATTCGGACCTTTGACTCCTGG"
        def cdr3Markup = new Cdr3Markup(7, 14, 19, 32)

        def pSegments = PSegmentSearcher.search(cdr3Markup, cdr3seq)

        assert pSegments.pSegmentEndV == -1
        assert pSegments.pSegmentStartJ == -1
    }

    @Test
    void okPTest() {
        //                    V1     D0   D1           J0
        //                    |    . |    |            |
        //             000000000011111111112222222222333333333344444
        //             012345678901234567890123456789012345678901234
        def cdr3seq = "TGTGCGAGCTCGCTTGGGGAAGACATTCGAAGCTTTGACTCCTGG"
        def cdr3Markup = new Cdr3Markup(7, 14, 19, 32)

        def pSegments = PSegmentSearcher.search(cdr3Markup, cdr3seq)

        assert pSegments.pSegmentEndV == 12
        assert pSegments.pSegmentStartJ == 26
    }

    @Test
    void twoConsMMPTest() {
        //                    V1     D0   D1           J0
        //                    |      |    |            |
        //             000000000011111111112222222222333333333344444
        //             012345678901234567890123456789012345678901234
        def cdr3seq = "TGTGCGAGGTGGCTTGGGGAAGGGAGAAAAAGCTTTGACTCCTGG"
        //                                   GGAGTCAAAG
        def cdr3Markup = new Cdr3Markup(7, 14, 19, 32)

        def pSegments = PSegmentSearcher.search(cdr3Markup, cdr3seq)

        assert pSegments.pSegmentStartJ == 28
    }

    @Test
    void singleBadSepMMPTest() {
        //                    V1     D0   D1           J0
        //                    |      |    |            |
        //             000000000011111111112222222222333333333344444
        //             012345678901234567890123456789012345678901234
        def cdr3seq = "TGTGCGAGGTGGCTTGGGGAAGGGTGCCAAAGCTTTGACTCCTGG"
        //                                   GGAGTCAAAG
        def cdr3Markup = new Cdr3Markup(7, 14, 19, 32)

        def pSegments = PSegmentSearcher.search(cdr3Markup, cdr3seq)

        assert pSegments.pSegmentStartJ == 27
    }
    
}
