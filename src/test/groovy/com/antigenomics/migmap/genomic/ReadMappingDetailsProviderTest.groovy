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

package com.antigenomics.migmap.genomic

import com.antigenomics.migmap.blast.BlastInstanceFactory
import com.antigenomics.migmap.io.Read
import com.antigenomics.migmap.mapping.ReadMappingDetailsProvider
import org.junit.AfterClass
import org.junit.Test

class ReadMappingDetailsProviderTest {
    private final BlastInstanceFactory factory

    ReadMappingDetailsProviderTest() {
        factory = new BlastInstanceFactory("data/", "human", ["IGH", "IGL"], true, false)
        factory.annotateV()
    }

    @Test
    void cdr23OverlapTest() {
        def instance = factory.create()

        def seq = "CTGAGTCAATCGCCCTCTGCCTCTGCCTCCCTGGGAGCCTCGGTCAAGCT" +
                "CCCCTGCACTGTCCAGTGGGCACGACTTCTACGCCATCGCATGCCATCAG" +
                "TAGCAGCCAGAGAAGGGGCCTCGCTTCTTGATGAAATCTAACAATGCTGG" +
                "CAGCCACAGTCAG"


        instance.put(new Read(null, seq))
        instance.put(null)

        def readMapping = instance.take()

        def details = ReadMappingDetailsProvider.getDetails(readMapping)

        def markup = readMapping.mapping.regionMarkup
        assert markup.cdr3Start < markup.cdr2End  // check for test case
        assert details.contignt.contains(seq)
    }

    @Test
    void fullLengthTest() {
        def instance = factory.create()

        def seq = "TAAGAGGGCAGTGGTATCAACGCAGAGTACGGATATTCTGAGGTCCGCTC" +
                "TCTTGGGGGGCTTTCTGAGAGTCGTGGATCTCATGTGCAAGAAAATGAAG" +
                "CACCTGTGGTTCTTCCTCCTGCTGGTGGCGGCTCCCAGATGGGTCCTGTC" +
                "CCAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCCTCGGAGA" +
                "CCCTGTCCCTCAGGTGCACTGTCTCTGGGGGCTCCATGACAAGAACTACT" +
                "GACTACTGGGGCTGGGTCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGAT" +
                "TGCAAGTGTCTCTTATAGTGGGAGCACCACCTACAACCCGTCCCGGAAGA" +
                "GTCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAA" +
                "CTGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTCTGTGCGAGGTG" +
                "GCTTGGGGAAGACATTCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGG" +
                "TCACCGTCTCTAA"

        instance.put(new Read(null, seq))
        instance.put(null)

        def readMapping = instance.take()

        def details = ReadMappingDetailsProvider.getDetails(readMapping)

        assert details.fr1nt == "CAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCCTCGGAGACCCTGTCCCTCAGGTGCACTGTCTCT"
        assert details.cdr1nt == "GGGGGCTCCATGACAAGAACTACTGACTAC"
        assert details.fr2nt == "TGGGGCTGGGTCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATTGCAAGT"
        assert details.cdr2nt == "GTCTCTTATAGTGGGAGCACC"
        assert details.fr3nt == "ACCTACAACCCGTCCCGGAAGAGTCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAACTGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTC"
        assert details.contignt == "CAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCCTCGGAGAC" +
                "CCTGTCCCTCAGGTGCACTGTCTCTGGGGGCTCCATGACAAGAACTACTG" +
                "ACTACTGGGGCTGGGTCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATT" +
                "GCAAGTGTCTCTTATAGTGGGAGCACCACCTACAACCCGTCCCGGAAGAG" +
                "TCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAAC" +
                "TGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTCTGTGCGAGGTGG" +
                "CTTGGGGAAGACATTCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGT" +
                "CACCGTCTCTAA"

        instance.close()
    }

    @Test
    void partialTest() {
        def instance = factory.create()

        def seq = "CTGGGTCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATTGCAAGTGTCT" +
                "CTTATAGTGGGAGCACCACCTACAACCCGTCCCGGAAGAGTCGAGTCACA" +
                "ATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAACTGAGGTCCAT" +
                "GACCGCCGCAGACACGGCTGTGTATTTCTGTGCGAGGTGGCTTGGGGAAG" +
                "ACATTCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTCT" +
                "AA"

        instance.put(new Read(null, seq))
        instance.put(null)

        def readMapping = instance.take()

        def details = ReadMappingDetailsProvider.getDetails(readMapping)

        assert details.fr1nt == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        assert details.cdr1nt == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        assert details.fr2nt == "NNNNNCTGGGTCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATTGCAAGT"
        assert details.cdr2nt == "GTCTCTTATAGTGGGAGCACC"
        assert details.fr3nt == "ACCTACAACCCGTCCCGGAAGAGTCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAACTGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTC"
        assert details.contignt == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                "NNNNNNNNNNCTGGGTCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATT" +
                "GCAAGTGTCTCTTATAGTGGGAGCACCACCTACAACCCGTCCCGGAAGAG" +
                "TCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAAC" +
                "TGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTCTGTGCGAGGTGG" +
                "CTTGGGGAAGACATTCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGT" +
                "CACCGTCTCTAA"

        instance.close()
    }

    @Test
    void partialTest2() {
        def instance = factory.create()

        def seq = "ATTATAGTGGGAGCACCACCTACAACCCGTCCCGGAAGAGTCGAGTCACA" +
                "ATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAACTGAGGTCCAT" +
                "GACCGCCGCAGACACGGCTGTGTATTTCTGTGCGAGGTGGCTTGGGGAAG" +
                "ACATTCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTCT" +
                "AA"

        instance.put(new Read(null, seq))
        instance.put(null)

        def readMapping = instance.take()

        def details = ReadMappingDetailsProvider.getDetails(readMapping)

        assert details.fr1nt == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        assert details.cdr1nt == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        assert details.fr2nt == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        assert details.cdr2nt == "NNNNATTATAGTGGGAGCACC"
        assert details.fr3nt == "ACCTACAACCCGTCCCGGAAGAGTCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAACTGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTC"
        assert details.contignt == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" +
                "NNNNNNNNNNATTATAGTGGGAGCACCACCTACAACCCGTCCCGGAAGAG" +
                "TCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAAC" +
                "TGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTCTGTGCGAGGTGG" +
                "CTTGGGGAAGACATTCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGT" +
                "CACCGTCTCTAA"

        assert details.contigaa == "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" +
                "XXXXYSGSTTYNPSRKSRVTISLDPSRNELSLELRSMTAADTAVYFCARW" +
                "LGEDIRTFDSWGQGTLVTVS"

        instance.close()
    }

    @AfterClass
    static void tearDown() {
        SegmentDatabase.clearTemporaryFiles()
    }
}
