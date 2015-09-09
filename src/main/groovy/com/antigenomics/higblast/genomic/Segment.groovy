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

class Segment {
    static final Segment DUMMY_J = new Segment(SegmentType.J, Util.MY_NA, "AAAAAAAAAAAAAAAAAAAAAAAAA", -1),
                         DUMMY_D = new Segment(SegmentType.D, Util.MY_NA, "AAAAAAAAAAAAAAAAAAAAAAAAA", -1)

    final String name, sequence, regexName
    final int referencePoint
    final SegmentType type

    Segment(SegmentType type, String name, String sequence, int referencePoint) {
        this.type = type
        this.name = name
        this.regexName = name.replace(".", "\\.").replace("*", "\\*")
        this.sequence = sequence
        this.referencePoint = referencePoint
    }

    String toFastaString() {
        ">$name\n$sequence"
    }

    int getFrame() {
        referencePoint % 3
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Segment segment = (Segment) o

        name == segment.name
    }

    @Override
    int hashCode() {
        name.hashCode()
    }


    @Override
    public String toString() {
        name
    }
}
