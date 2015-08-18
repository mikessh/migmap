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

import com.antigenomics.higblast.Util

class JSegment extends Segment {
    static final JSegment DUMMY = new JSegment(Util.MY_NA, "AAAAAAAAAAAAAAAAAAAAAAAAA", -1)

    JSegment(String name, String sequence, int referencePoint) {
        super(SegmentType.J,name, sequence, referencePoint)
    }
}
