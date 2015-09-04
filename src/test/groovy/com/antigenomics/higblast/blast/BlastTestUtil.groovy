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
import com.antigenomics.higblast.io.Read

class BlastTestUtil {
    final static SegmentDatabase segmentDatabase = new SegmentDatabase("data/", "human", ["IGH"])
    final static BlastParser parser = new BlastParser(segmentDatabase)

    static Read toRead(String seq) {
        new Read("@test", seq, "I" * seq.length())
    }

    static String toChunk(Read read) {
        def factory = new BlastInstanceFactory("data/", "human", ["IGH"], true, false)
        def instance = factory.create()

        instance.put(read)
        instance.put(null)

        instance.nextChunk()
    }
}
