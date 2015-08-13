package com.antigenomics.higblast.shm
/**
 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */
class Hypermutation {
    final int pos, posInRead
    final String ntFrom, ntTo, aaFrom, aaTo
    final String region

    Hypermutation(int pos, int posInRead,
                  char ntFrom, char ntTo, char aaFrom, char aaTo,
                  String region) {
        this.pos = pos
        this.posInRead = posInRead
        this.ntFrom = ntFrom
        this.ntTo = ntTo
        this.aaFrom = aaFrom
        this.aaTo = aaTo
        this.region = region
    }

    boolean equals(o) {
        if (this.is(o)) return true

        Hypermutation that = (Hypermutation) o

        if (ntFrom != that.ntFrom) return false
        if (ntTo != that.ntTo) return false
        if (pos != that.pos) return false

        true
    }

    int hashCode() {
        int result
        result = pos
        result = 31 * result + (int) ntFrom
        31 * result + (int) ntTo
    }

    String toString() {
        "$pos:$ntFrom>$ntTo,${(int) (pos / 3)}:$aaFrom>$aaTo,$region"
    }
}
