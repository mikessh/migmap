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

import com.antigenomics.higblast.io.Read
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

        def extractedRead = instance.getRead(chunk)
        assert extractedRead.equals(read)

        assert instance.nextChunk() == null
        assert chunk.contains("# V-(D)-J rearrangement summary")
        assert chunk.contains("# V-(D)-J junction details")
        assert chunk.contains("# Alignment summary between query and top germline V gene hit")
        assert chunk.contains("# Hit table ")

        println chunk

        assert instance.proc.exitValue() == 0

        assert instance.parse(chunk).cdr3nt == "TGTGCGAGGTGGCTTGGGGAAGACATTCGGACCTTTGACTCCTGG"

        factory.segmentDatabase.clearBlastDb()
    }

    @Test
    void multiQueryTest() {
        int nQueries = 100
        def factory = new BlastInstanceFactory("data/", "human", ["IGH"], true, false)
        def instance = factory.create()

        def readsIds = new HashSet<String>()

        (0..<nQueries).each {
            def header = "@$it"
            readsIds.add(header)
            instance.put(new Read(header, seq, qual))
        }

        instance.put(null)

        def extractedIds = new HashSet<String>((0..<nQueries).collect {
            def chunk = instance.nextChunk()
            assert instance.parse(chunk).cdr3nt == "TGTGCGAGGTGGCTTGGGGAAGACATTCGGACCTTTGACTCCTGG"
            BlastInstance.getRead(chunk).header
        })

        def intersection = extractedIds.intersect(readsIds)

        assert extractedIds.size() == nQueries
        assert intersection.size() == nQueries

        assert instance.proc.exitValue() == 0

        factory.segmentDatabase.clearBlastDb()
    }
}
