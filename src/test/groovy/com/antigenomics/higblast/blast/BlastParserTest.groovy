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

package com.antigenomics.higblast.blast

import com.antigenomics.higblast.genomic.SegmentDatabase
import org.junit.AfterClass
import org.junit.Test

class BlastParserTest {
    @Test
    void exactTest() {
        def segmentDatabase = new SegmentDatabase("data/", "human", ["IGH"], true, false)
        def parser = new BlastParser(segmentDatabase)

        def chunk = "# IGBLASTN 2.2.29+\n" +
                "# Query: @MIG UMI:GGATATGCCGCTC:8|TAAGAGGGCAGTGGTATCAACGCAGAGTACGGATATTCTGAGGTCCGCTCTCTTGGGGGGCTTTCTGAGAGTCGTGGATCTCATGTGCAAGAAAATGAAGCACCTGTGGTTCTTCCTCCTGCTGGTGGCGGCTCCCAGATGGGTCCTGTCCCAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCCTCGGAGACCCTGTCCCTCAGGTGCACTGTCTCTGGGGGCTCCATGACAAGAACTACTGACTACTGGGGCTGGGTCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATTGCAAGTGTCTCTTATAGTGGGAGCACCACCTACAACCCGTCCCGGAAGAGTCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAACTGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTCTGTGCGAGGTGGCTTGGGGAAGACATTCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTCTAA|#'.IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIBIIIIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIB;IIIIIIIIIIII.''\n" +
                "# Database: data//database-66ede6f0-a68a-4931-bd32-5cabd404d955/v data//database-66ede6f0-a68a-4931-bd32-5cabd404d955/d data//database-66ede6f0-a68a-4931-bd32-5cabd404d955/j\n" +
                "# Domain classification requested: imgt\n" +
                "\n" +
                "# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.\n" +
                "IGHV4-39*01\tIGHD6-19*01\tIGHJ4*02\tVH\tNo\tN/A\tN/A\t+\n" +
                "\n" +
                "# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself\n" +
                "GCGAG\tN/A\tGTGGCT\tTGGGGAAGACATTCGGAC\tCTTTG\t\n" +
                "\n" +
                "# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)\n" +
                "FR1-IMGT\t152\t226\t75\t72\t3\t0\t96\n" +
                "CDR1-IMGT\t227\t256\t30\t22\t8\t0\t73.3\n" +
                "FR2-IMGT\t257\t307\t51\t47\t4\t0\t92.2\n" +
                "CDR2-IMGT\t308\t328\t21\t19\t2\t0\t90.5\n" +
                "FR3-IMGT\t329\t442\t114\t96\t18\t0\t84.2\n" +
                "CDR3-IMGT (germline)\t443\t447\t5\t5\t0\t0\t100\n" +
                "Total\tN/A\tN/A\t296\t261\t35\t0\t88.2\n" +
                "\n" +
                "# Hit table (the first field indicates the chain type of the hit)\n" +
                "# Fields: query id, q. start, query seq, s. start, subject seq\n" +
                "# 3 hits found\n" +
                "V\t@MIG\t152\tCAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCCTCGGAGACCCTGTCCCTCAGGTGCACTGTCTCTGGGGGCTCCATGACAAGAACTACTGACTACTGGGGCTGGGTCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATTGCAAGTGTCTCTTATAGTGGGAGCACCACCTACAACCCGTCCCGGAAGAGTCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAACTGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTCTGTGCGAG\t1\tCAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGCAGTAGTAGTTACTACTGGGGCTGGATCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGAGTATCTATTATAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCCGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCAGACACGGCTGTGTATTACTGTGCGAG\n" +
                "D\t@MIG\t448\tGTGGCT\t11\tGTGGCT\n" +
                "J\t@MIG\t472\tCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTC\t5\tCTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTC"

        def mapping = parser.parse(chunk)

        assert mapping.vSegments[0].name == "IGHV4-39*01"
        assert mapping.dSegments[0].name == "IGHD6-19*01"
        assert mapping.jSegments[0].name == "IGHJ4*02"

        assert mapping.cdr3Markup.vEnd == 8
        assert mapping.cdr3Markup.dStart == 8
        assert mapping.cdr3Markup.dEnd == 14
        assert mapping.cdr3Markup.jStart == 32

        assert mapping.regionMarkup.cdr1Start == 226
        assert mapping.regionMarkup.cdr1End == 256
        assert mapping.regionMarkup.cdr2Start == 307
        assert mapping.regionMarkup.cdr2End == 328
    }

    @Test
    void mutationCoordinatesTest() {
        // def seq = "GACACGGCTGTGTATTTCTGTGCGAGATATAGTGACTACGATTACCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTC"
        //            ................A.........--------------------------------------------------------------- IGHV3-33*01
        //            --------------------------........T..........-------------------------------------------- IGHD5/OR15-5a*01
        //            --------------------------------------------------........A.............................. IGHJ4*02
        //            78901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
        //            66677777777778888888888999999999900000000001111111111222222222233333333334444444444555555
        //            22222222222222222222222222222222233333333333333333333333333333333333333333333333333333333

        // def factory = new BlastInstanceFactory("data/", "human", ["IGH"], true, false)
        // def instance = factory.create()

        // def read = new Read(
        //        "@",
        //        seq,
        //        "I" * seq.length()
        // )

        // instance.put(read)
        // instance.put(null)
        // def readMapping = instance.take()

        def chunk = "# IGBLASTN 2.2.29+\n" +
                "# Query: @|GACACGGCTGTGTATTTCTGTGCGAGATATAGTGACTACGATTACCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTC|IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n" +
                "# Database: data//database-a52bccf8-da3b-473f-bacc-f9d34e532516/v data//database-a52bccf8-da3b-473f-bacc-f9d34e532516/d data//database-a52bccf8-da3b-473f-bacc-f9d34e532516/j\n" +
                "# Domain classification requested: imgt\n" +
                "\n" +
                "# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.\n" +
                "IGHV3-33*01\tIGHD5/OR15-5a*01\tIGHJ4*02\tVH\tNo\tN/A\tN/A\t+\n" +
                "\n" +
                "# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself\n" +
                "TGCGA\t(GA)\tTATAGTGACTACGATTAC\tCGGAC\tCTTTG\t\n" +
                "\n" +
                "# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)\n" +
                "FR3-IMGT\t1\t21\t21\t20\t1\t0\t95.2\n" +
                "CDR3-IMGT (germline)\t22\t27\t6\t6\t0\t0\t100\n" +
                "Total\tN/A\tN/A\t27\t26\t1\t0\t96.3\n" +
                "\n" +
                "# Hit table (the first field indicates the chain type of the hit)\n" +
                "# Fields: query id, q. start, query seq, s. start, subject seq\n" +
                "# 3 hits found\n" +
                "V\t@|GACACGGCTGTGTATTTCTGTGCGAGATATAGTGACTACGATTACCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTC|IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\t1\tGACACGGCTGTGTATTTCTGTGCGAGA\t268\tGACACGGCTGTGTATTACTGTGCGAGA\n" +
                "D\t@|GACACGGCTGTGTATTTCTGTGCGAGATATAGTGACTACGATTACCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTC|IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\t26\tGATATAGTGACTACGATTAC\t4\tGATATAGTGTCTACGATTAC\n" +
                "J\t@|GACACGGCTGTGTATTTCTGTGCGAGATATAGTGACTACGATTACCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTC|IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\t51\tCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTC\t5\tCTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTC\n" +
                "# BLAST processed 1 queries"

        def segmentDatabase = new SegmentDatabase("data/", "human", ["IGH"], true, false)
        def parser = new BlastParser(segmentDatabase)

        def mapping = parser.parse(chunk)
        def mutations = mapping.mutations

        assert mutations.size() == 3
        assert mutations[0].pos == 283
        assert mutations[1].pos == 301
        assert mutations[2].pos == 325
    }

    @AfterClass
    static void tearDown() {
        SegmentDatabase.clearTemporaryFiles()
    }
}
