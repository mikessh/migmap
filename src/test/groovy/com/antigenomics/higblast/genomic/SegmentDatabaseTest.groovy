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

package com.antigenomics.higblast.genomic

import com.antigenomics.higblast.blast.BlastInstanceFactory
import org.junit.AfterClass
import org.junit.Test

class SegmentDatabaseTest {
    @Test
    void massiveLoadTest() {
        ["human", "mouse", "rat", "rabbit", "rhesus_monkey"].each {
            def segmentDatabase = new SegmentDatabase("data/", it, ["TRA", "TRB", "TRG", "TRD",
                                                                    "IGH", "IGK", "IGL"])

            assert !segmentDatabase.segments.isEmpty()
            assert segmentDatabase.vSegments > 0
            assert segmentDatabase.dSegments > 0
            assert segmentDatabase.jSegments > 0
        }
    }

    @Test
    void loadTest() {
        def segmentDatabase = new SegmentDatabase("data/", "human", ["IGH"])

        assert !segmentDatabase.segments.isEmpty()

        def testSegments = ["IGHV4-59*09",
                            "IGHV3-23*01",
                            "IGHV3/OR16-8*01",
                            "IGHV3-73*02",
                            "IGHV4-39*06",
                            "IGHV3-30-3*03",
                            "IGHV1/OR15-2*01",
                            "IGHV4/OR15-8*02",
                            "IGHV4-28*05",
                            "IGHV1/OR15-3*02"]

        assert testSegments.each { segmentDatabase.segments.keySet().contains(it) }
    }

    @Test
    void makeDbTest() {
        def segmentDatabase = new SegmentDatabase("data/", "human", ["TRA"])
        segmentDatabase.makeBlastDb()

        assert new File(segmentDatabase.databaseTempPath).exists()
        assert new File("$segmentDatabase.databaseTempPath/v.nhr").exists()
        assert new File("$segmentDatabase.databaseTempPath/d.nhr").exists()
        assert new File("$segmentDatabase.databaseTempPath/j.nhr").exists()
    }

    @Test
    void annotationTestIG() {
        def factory = new BlastInstanceFactory("data/", "human", ["IGH", "IGK", "IGL"], true, false)
        factory.annotateV()

        assert factory.segmentDatabase.annotatedV == factory.segmentDatabase.vSegments
    }

    @Test
    void annotationTestTR() {
        def factory = new BlastInstanceFactory("data/", "human", ["TRA", "TRB", "TRG", "TRD"], true, false)
        factory.annotateV()

        assert factory.segmentDatabase.annotatedV + 2 == factory.segmentDatabase.vSegments // known 2 unannotateable cases
    }

    @AfterClass
    static void tearDown() {
        SegmentDatabase.clearTemporaryFiles()
    }
}
