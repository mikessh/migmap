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

package com.antigenomics.higblast.mutation

import com.antigenomics.higblast.blast.Alignment
import com.antigenomics.higblast.genomic.SegmentDatabase
import com.antigenomics.higblast.mapping.RegionMarkup
import org.junit.AfterClass
import org.junit.Test

import static com.antigenomics.higblast.mutation.SubRegion.*

class MutationExtractorTest {
    @Test
    void test() {
        def alignment = new Alignment(0,
                //0000000001  11111111122222222223333
                //1234567890  12345678901234567890123
                "AGATCGATCGA--CTGCTACGACTGCATGACTCAAT", 0,
                //0000000001111111111222222   2222333
                //1234567890123456789012345   6789012
                "AAATCGATCGAAACTGCTACGACTGC---ACTCAGT")

        def posInReadExp = [1, 11, 24, 32]
        def expectedCodes = ["S1:A>G", "D11:AA", "I26:ATG", "S31:G>A"]
        def mutations = MutationExtractor.extract(alignment).mutations

        println mutations

        assert expectedCodes.size() == mutations.size()

        mutations.eachWithIndex { it, i ->
            assert it.toString() == expectedCodes[i]
            assert it.posInRead == posInReadExp[i]
        }
    }

    @Test
    void markupTest() {
        def query = "CAGCTGGAGTTGGTACAGTCTGGGGCTGAGGAGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATCCACATTCAGC" +
                "GGCCACTTTATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGGTGGATCAACTCTTACAGTGGTGCCACAAAGTAT" +
                "GCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATGACCACAATCTACATGGAGCTGAGCGGACTCACATCTGACGAC" +
                "ACGGCCGTGTATTTTTGTACCAGA"

        def segmentDatabase = new SegmentDatabase("data/", "human", ["IGH"])

        def segment = segmentDatabase.segments["IGHV1-2*02"]
        def regionMarkup = new RegionMarkup(75, 99, 150, 174, 287, 288)
        def alignment = new Alignment(0, query, 0, segment.sequence.substring(0, query.length()))
        def mutations = new MutationExtractor(segment, alignment, regionMarkup).mutations

        def expectedSubRegions = [FR1, FR1, FR1, FR1, FR1,
                                  CDR1, CDR1, CDR1, CDR1, CDR1,
                                  FR2,
                                  CDR2, CDR2, CDR2,
                                  FR3, FR3, FR3, FR3, FR3, FR3, FR3, FR3, FR3, FR3, FR3,
                                  CDR3, CDR3]

        mutations.eachWithIndex { it, i ->
            assert it.subRegion == expectedSubRegions[i]
        }
    }

    @AfterClass
    static void tearDown() {
        SegmentDatabase.clearTemporaryFiles()
    }
}
