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

package com.antigenomics.migmap

import com.antigenomics.migmap.blast.BlastInstanceFactory
import com.antigenomics.migmap.genomic.SegmentDatabase
import com.antigenomics.migmap.io.*
import com.antigenomics.migmap.mapping.ReadMapping
import com.antigenomics.migmap.mapping.ReadMappingDetailsProvider
import com.antigenomics.migmap.mapping.ReadMappingFilter
import org.junit.AfterClass
import org.junit.Test

import java.util.concurrent.atomic.AtomicInteger

class PipelineTest {
    private final BlastInstanceFactory factory = new BlastInstanceFactory("data/", "human", ["IGH"], true, false)

    PipelineTest() {
        factory.annotateV()
    }

    private static class NullOutputStream extends OutputStream {
        static final NullOutputStream INSTANCE = new NullOutputStream()

        @Override
        void write(int b) throws IOException {
        }
    }

    @Test
    void sampleTest() {
        def reader = new FastqReader("sample.fastq.gz", true)
        def filter = new ReadMappingFilter()

        def pipeline = new Pipeline(reader, factory,
                new InputPortMerge(
                        new ReadMappingOutput(new PlainTextOutput(NullOutputStream.INSTANCE),
                                new ReadMappingDetailsProvider()),
                        new ClonotypeOutput(new PlainTextOutput(NullOutputStream.INSTANCE))),
                filter)

        pipeline.run()

        assert pipeline.inputCount == 1000
        assert pipeline.readMappingFilter.mappedRatio >= 0.95
        assert pipeline.readMappingFilter.noCdr3Ratio <= 0.1
        assert pipeline.readMappingFilter.incompleteRatio <= 0.1
        assert pipeline.readMappingFilter.nonCanonicalRatio <= 0.1
    }

    @Test
    void fastaTest() {
        def reader = new FastaReader("sample.fasta.gz", true)

        def pipeline = new Pipeline(reader, factory,
                DummyInputPort.INSTANCE,
                ReadMappingFilter.createDummy())

        pipeline.run()

        assert pipeline.readMappingFilter.total == 100
        assert pipeline.readMappingFilter.goodRatio == 1.0
    }

    @Test
    void dAssignmentTest() {
        def reader = new FastqReader("ambiguous_d.fastq.gz", true)

        def badDCount = new AtomicInteger(),
            wrongMappingCount = new AtomicInteger()

        def pipeline = new Pipeline(reader, factory,
                new InputPort<ReadMapping>() {
                    @Override
                    void put(ReadMapping obj) {
                        if (obj.mapped && obj.cdr3nt == "TGTGCGAGCGATCGGAACGGTATGGACGTCTGG") {
                            if (!obj.mapping.dSegment.name == "IGHD1-1*01")
                                badDCount.incrementAndGet()
                        } else {
                            wrongMappingCount.incrementAndGet()
                        }
                    }

                    @Override
                    void close() {

                    }
                },
                ReadMappingFilter.createDummy())

        pipeline.run()

        assert badDCount.get() == 0
        assert wrongMappingCount.get() == 0
    }

    @Test
    void badDataTest() {
        def reader = new FastqReader("bad_sample.fastq.gz", true)

        def filter = new ReadMappingFilter((byte) 20, true, true, true, true)

        def pipeline = new Pipeline(reader, factory,
                new ReadMappingOutput(new PlainTextOutput(NullOutputStream.INSTANCE),
                        new ReadMappingDetailsProvider()),
                filter)

        pipeline.run()

        assert filter.passed == filter.good
        assert pipeline.readMappingFilter.mappedRatio >= 0.90
        assert 1 - pipeline.readMappingFilter.noCdr3Ratio >= 0.03
    }

    @AfterClass
    static void tearDown() {
        SegmentDatabase.clearTemporaryFiles()
    }
}
