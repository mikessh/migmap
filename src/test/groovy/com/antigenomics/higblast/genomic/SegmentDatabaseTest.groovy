/*
 * Copyright 2013-{year} Mikhail Shugay (mikhail.shugay@gmail.com)
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

package com.antigenomics.higblast.genomic

import org.junit.Test

class SegmentDatabaseTest {
    @Test
    void massiveLoadTest() {
        SegmentDatabase.SPECIES_ALIAS.keySet().each {
            println "Loading data for $it"
            
            def segmentDatabase = new SegmentDatabase("data/", it, ["TRA", "TRB", "TRG", "TRD",
                                                          "IGH", "IGK", "IGL"], true, false)

            assert !segmentDatabase.segments.isEmpty()
        }
    }

    @Test
    void loadTest() {
        def segmentDatabase = new SegmentDatabase("data/", "human", ["IGH"], true, false)

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
        def segmentDatabase = new SegmentDatabase("data/", "human", ["TRA"], true, false)
        segmentDatabase.makeBlastDb()

        assert new File(segmentDatabase.databaseTempPath).exists()
        assert new File("$segmentDatabase.databaseTempPath/v.nhr").exists()
        assert new File("$segmentDatabase.databaseTempPath/d.nhr").exists()

        segmentDatabase.clearBlastDb()
    }

    @Test
    void markupTest() {
        def segmentDatabase = new SegmentDatabase("data/", "human", ["IGH"], true, false)
        def segment = segmentDatabase.segments["IGHV1-18*01"]

        assert (segment as VSegment).cdr1start == 75
        assert (segment as VSegment).cdr1end == 99
        assert (segment as VSegment).cdr2start == 150
        assert (segment as VSegment).cdr2end == 174
    }
}
