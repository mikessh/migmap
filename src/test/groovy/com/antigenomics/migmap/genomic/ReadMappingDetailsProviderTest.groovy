/*
 * Copyright 2014-2015 Mikhail Shugay
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.migmap.genomic

import com.antigenomics.migmap.blast.BlastInstanceFactory
import com.antigenomics.migmap.io.Read
import com.antigenomics.migmap.mapping.ReadMappingDetailsProvider
import org.junit.AfterClass
import org.junit.Test

class ReadMappingDetailsProviderTest {
    private final BlastInstanceFactory factory, factoryMm

    ReadMappingDetailsProviderTest() {
        factory = new BlastInstanceFactory("data/", "human", ["IGH", "IGL"], true, false)
        factory.annotateV()
        factoryMm = new BlastInstanceFactory("data/", "mouse", ["IGH"], true, false)
        factoryMm.annotateV()
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
    void rcTest() {
        def instance = factoryMm.create()

        def seq = "CAGATCCTCTTCTGAGATGAGTTTTTGTTCGCTACCGCCACCCTCGAGTG" +
                "AGGAGACGGTGACCATTGTCCCTTGGCCCCAGATATCAAAAGCACTCCTA" +
                "GTTCCAGTTAAACCTCCTCTTGTACAATAATACACAGCCGTGTCCTCGGG" +
                "AGTCACAGAGTTCAGCTGCAGGGAGAACTGGTTCTTGGATGTGTCTGGGT" +
                "TGATGCTTATTCGACTTTTCACAGACACTGCATATTCATTATACCACTTG" +
                "GACCTGTAGTATATCCTTCCCAGCCACTCAAGGCCTCTCGATGGGGACTG" +
                "CCTGATCCAGTTCCAAGCAGCACTGTTGCTAGAGACACTGTCCCCGGAGA" +
                "TGGCACAGGTGAGTAAGAGGGTCTGCGAGGGCTTCACCAGTCCTGGACCT" +
                "GACTGCTGCAGCTGTACCTGGGCCATGGCCGGCTGGGCCGCGAGTAATAA" +
                "CAATCCAGCGGC"

        instance.put(new Read(null, seq))
        instance.put(null)

        def readMapping = instance.take()

        def details = ReadMappingDetailsProvider.getDetails(readMapping)

        assert details.cdr1aa == "GDSVSSNSAA"
        assert details.cdr2aa == "IYYRSKWYN"

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
