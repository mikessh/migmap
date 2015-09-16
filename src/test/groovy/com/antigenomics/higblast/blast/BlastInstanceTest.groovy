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

package com.antigenomics.higblast.blast

import com.antigenomics.higblast.genomic.SegmentDatabase
import com.antigenomics.higblast.io.Read
import org.junit.AfterClass
import org.junit.Test

class BlastInstanceTest {
    final String seq = "TAAGAGGGCAGTGGTATCAACGCAGAGTACGGATATTCTGAGGTCCGCTC" +
            "TCTTGGGGGGCTTTCTGAGAGTCGTGGATCTCATGTGCAAGAAAATGAAG" +
            "CACCTGTGGTTCTTCCTCCTGCTGGTGGCGGCTCCCAGATGGGTCCTGTC" +
            "CCAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCCTCGGAGA" +
            "CCCTGTCCCTCAGGTGCACTGTCTCTGGGGGCTCCATGACAAGAACTACT" +
            "GACTACTGGGGCTGGGTCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGAT" +
            "TGCAAGTGTCTCTTATAGTGGGAGCACCACCTACAACCCGTCCCGGAAGA" +
            "GTCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAA" +
            "CTGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTCTGTGCGAGGTG" +
            "GCTTGGGGAAGACATTCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGG" +
            "TCACCGTCTCTAA",
                 qual = "#'.IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" +
                         "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" +
                         "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" +
                         "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIB" +
                         "IIIIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIII" +
                         "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" +
                         "IIIIIIIIIIIBIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" +
                         "IIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" +
                         "IIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" +
                         "IIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIB;II" +
                         "IIIIIIIIII.''"

    @Test
    void singleQueryTest() {
        def factory = new BlastInstanceFactory("data/", "human", ["IGH"], true, false)
        def instance = factory.create()

        def read = new Read(
                "@MIG UMI:GGATATGCCGCTC:8",
                seq,
                qual
        )

        instance.put(read)

        instance.put(null)

        //instance.proc.in.eachLine { println it }

        String chunk = instance.nextChunk()

        //def extractedRead = instance.getRead(chunk)
        //assert extractedRead.equals(read)
        assert instance.getHeader(chunk).equals(read.header)

        assert instance.nextChunk() == null
        assert chunk.contains("# V-(D)-J rearrangement summary")
        assert chunk.contains("# V-(D)-J junction details")
        assert chunk.contains("# Alignment summary between query and top germline V gene hit")
        assert chunk.contains("# Hit table ")

        println chunk

        assert instance.proc.exitValue() == 0
    }

    @Test
    void multiQueryTest() {
        int nQueries = 100
        def factory = new BlastInstanceFactory("data/", "human", ["IGH"], true, false)
        def instance = factory.create()

        def readsIds = new HashSet<String>()

        (0..<nQueries).each {
            def header = "@$it".toString()
            readsIds.add(header)
            instance.put(new Read(header, seq, qual))
        }

        instance.put(null)

        def extractedIds = new HashSet<String>((0..<nQueries).collect {
            def mapping = instance.take()
            assert mapping.cdr3nt == "TGTGCGAGGTGGCTTGGGGAAGACATTCGGACCTTTGACTCCTGG"
            mapping.read.header
        })

        def intersection = extractedIds.intersect(readsIds)

        assert extractedIds.size() == nQueries
        assert intersection.size() == nQueries

        assert instance.proc.exitValue() == 0
    }

    @AfterClass
    static void tearDown() {
        SegmentDatabase.clearTemporaryFiles()
    }
}
