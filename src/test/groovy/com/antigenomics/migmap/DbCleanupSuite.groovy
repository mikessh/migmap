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

package com.antigenomics.migmap

import com.antigenomics.migmap.blast.BlastInstanceFactoryTest
import com.antigenomics.migmap.blast.BlastParserTest
import com.antigenomics.migmap.blast.DSegmentTest
import com.antigenomics.migmap.blast.PSegmentTest
import com.antigenomics.migmap.blast.RefPointSearcherTest
import com.antigenomics.migmap.genomic.ReadMappingDetailsProviderTest
import com.antigenomics.migmap.genomic.SegmentDatabase
import com.antigenomics.migmap.genomic.SegmentDatabaseTest
import org.junit.AfterClass
import org.junit.runner.RunWith
import org.junit.runners.Suite

@RunWith(Suite.class)
@Suite.SuiteClasses([PipelineTest.class,
        BlastInstanceFactoryTest.class,
        BlastParserTest.class,
        DSegmentTest.class,
        PSegmentTest.class,
        RefPointSearcherTest.class,
        ReadMappingDetailsProviderTest.class,
        SegmentDatabaseTest.class])
class DbCleanupSuite {
    @AfterClass
    static void tearDown() {
        SegmentDatabase.clearTemporaryFiles()
    }
}