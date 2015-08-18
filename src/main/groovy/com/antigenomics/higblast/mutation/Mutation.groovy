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

package com.antigenomics.higblast.mutation

import com.antigenomics.higblast.genomic.Segment

class Mutation {
    final int startInRead, endInRead
    final String ntFrom, ntTo

    int start, end

    String aaFrom, aaTo // todo

    final MutationType type
    Segment region
    SubRegion subRegion

    Mutation(MutationType type,
             int start, int end, int startInRead, int endInRead,
             String ntFrom, String ntTo) {
        this.type = type
        this.end = end
        this.startInRead = startInRead
        this.endInRead = endInRead
        this.ntFrom = ntFrom
        this.ntTo = ntTo
        this.start = start
    }

    int getPos() {
        type == MutationType.Insertion ? end : start
    }

    int getPosInRead() {
        type == MutationType.Deletion ? endInRead : startInRead // deletion in ref == insertion in query
    }

    @Override
    public String toString() {
        // insertion occur before "position", so that seq[position, ...] is shifted forward
        type.shortName + pos + ":" +
                (type == MutationType.Insertion ? "" : ntFrom) +
                (type == MutationType.Substitution ? ">" : "") +
                (type == MutationType.Deletion ? "" : ntTo)
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Mutation mutation = (Mutation) o

        start == mutation.start &&
                ntFrom == mutation.ntFrom &&
                ntTo == mutation.ntTo &&
                region == mutation.region
    }

    @Override
    int hashCode() {
        int result
        result = ntFrom.hashCode()
        result = 31 * result + ntTo.hashCode()
        result = 31 * result + start
        31 * result + region.hashCode()
    }
}
