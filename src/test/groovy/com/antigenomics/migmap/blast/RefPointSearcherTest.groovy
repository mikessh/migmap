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

import com.antigenomics.migmap.genomic.SegmentDatabase
import org.junit.AfterClass
import org.junit.Test

import static com.antigenomics.migmap.blast.BlastTestUtil.*

class RefPointSearcherTest {
    final BlastInstance blastInstance

    RefPointSearcherTest() {
        blastInstance = new BlastInstanceFactory("data/", "human", ["IGH"], true, false).create()
    }

    @Test
    void cdr3EndTest() {
        def seq = "CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTAAAGCCTTCGGAGACCCTGTCCCTCACCTGCGGTGTCTCTGGTTACTCCATAAGCAGTGGTTACTACTGGGCCTGGATCCGGCAGCCCCCAGGGAAGGGGCTGGAGTGGGTTGCGACTATCTATCATGATGGAAGATCCTACTACAACCCGTCCCTGGAAAGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAGGCTGACTTCTGTGACCGCCGCAGACACGGCCGTATATTTCTGTGCGAGGGATCAGGGGCCACGAGACCACCCCAGTTCATCGTTTTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA"
        def read = toRead(seq)
        def chunk = toChunk(read)

        def mapping = parser.parse(chunk)
        assert mapping.hasCdr3
        assert mapping.complete
        assert mapping.regionMarkup.cdr3End == 333

        def readMapping = blastInstance.createReadMapping(mapping, read)
        assert readMapping.canonical
        assert readMapping.cdr3nt == "TGTGCGAGGGATCAGGGGCCACGAGACCACCCCAGTTCATCGTTTTGG"
    }

    @Test
    void cdr3StartTest() {
        def seq = "CAGCTGCGGCTGCAGGAGTCGGGCCCAGGACAGGTGAAGCCTTCGGAGACCCTGTCCCTCACTTGCACGGTCTCTGCTGGCGCCATCACCAGTGATCATTACTACTGGGGCTGGGTCCGCCAGCGCCCAGGGAAGGGACTAGAGTGGATTGGGAGTGTCCATAATAGTGGGAGCACCTCCTACAACCCGTCCCTTCAGAGTCGAATCACCATGTCCATAGACATGTCGAAGAACCACTTCTCCCTGCGGCTGACTTCCGTGACCGCCGCAGACACGGGTGTTTATTACTGCACCACCTATGGTTATTATTACTATTACGGCATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCA"
        def read = toRead(seq)
        def chunk = toChunk(read)

        def mapping = parser.parse(chunk)

        def readMapping = blastInstance.createReadMapping(mapping, read)
        assert readMapping.cdr3nt == "TGCACCACCTATGGTTATTATTACTATTACGGCATGGACGTCTGG"
    }

    @Test
    void cdr3StartTest2() {
        def seq = "CAGGTGCAGCTGATGGAGTCTGGGGGAGGCGTTGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGAAGCCTCTGGATTCACCTTCAGTCGTTCTTCTATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGTTGGCACTTATTTCACGTGATGGGAACTATGAGAGTTACGCAGACTCCGTGAAGGGCCGATTCTCGATCTCCAGAGACAACTCCAAGAACACTCTCTTTCTGCAAGTGGACGGCCTGAGAGCTGAGGACACGGCTGTGTATTATTGTCTTGGTGACAACGGCCATTGGGGCCAGGGAACCCTGCTCAGCGTCTCCTCA"
        def read = toRead(seq)
        def chunk = toChunk(read)

        def mapping = parser.parse(chunk)

        def readMapping = blastInstance.createReadMapping(mapping, read)
        assert readMapping.cdr3nt == "TGTCTTGGTGACAACGGCCATTGG"
    }

    @Test
    void borderlineCase() {
        def chunk = "# IGBLASTN 2.2.29+\n" +
                "# Query: @M01691:112:000000000-AFF3F:1:1106:14463:20458 2:N:0:1 UMI:ATGTGTACCGGG:GGGGGGGGGGGG\n" +
                "# Database: /Users/mikesh/Programming/migmap/data/database-f1894158-8d79-41c0-8288-0bd8217f7c4f/v /Users/mikesh/Programming/migmap/data/database-f1894158-8d79-41c0-8288-0bd8217f7c4f/d /Users/mikesh/Programming/migmap/data/database-f1894158-8d79-41c0-8288-0bd8217f7c4f/j\n" +
                "# Domain classification requested: imgt\n" +
                "\n" +
                "# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.\n" +
                "IGHV3-21*01\tIGHD5-5*01\tIGHJ6*04\tVH\tYes\tN/A\tNo\t+\n" +
                "\n" +
                "# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself\n" +
                "GAGAG\tTCCCAGAGGC\tGGATACAGCTATGGTT\tTTGAG\tTTACT\t\n" +
                "\n" +
                "# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)\n" +
                "FR1-IMGT\t93\t167\t75\t75\t0\t0\t100\n" +
                "CDR1-IMGT\t168\t191\t24\t24\t0\t0\t100\n" +
                "FR2-IMGT\t192\t242\t51\t51\t0\t0\t100\n" +
                "CDR2-IMGT\t243\t266\t24\t24\t0\t0\t100\n" +
                "FR3-IMGT\t267\t376\t114\t110\t0\t4\t96.5\n" +
                "CDR3-IMGT (germline)\t377\t383\t7\t7\t0\t0\t100\n" +
                "Total\tN/A\tN/A\t295\t291\t0\t4\t98.6\n" +
                "\n" +
                "# Hit table (the first field indicates the chain type of the hit)\n" +
                "# Fields: subject id, q. start, query seq, s. start, subject seq\n" +
                "# 3 hits found\n" +
                "V\tIGHV3-21*01\t93\tGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCCTGGTCAAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCATCCATTAGTAGTAGTAGTAGTTACATATACTAC----ACTCAGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAG\t1\tGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCCTGGTCAAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCATCCATTAGTAGTAGTAGTAGTTACATATACTACGCAGACTCAGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAG\n" +
                "D\tIGHD5-5*01\t394\tGGATACAGCTATGGTT\t3\tGGATACAGCTATGGTT\n" +
                "J\tIGHJ6*04\t415\tTTACTACTACTACTACGGTATGGACGTCTGG\t2\tTTACTACTACTACTACGGTATGGACGTCTGG\n" +
                "# BLAST processed 1 queries"

        def mapping = blastInstance.parser.parse(chunk)

        assert mapping.complete
    }

    @AfterClass
    static void tearDown() {
        SegmentDatabase.clearTemporaryFiles()
    }
}
