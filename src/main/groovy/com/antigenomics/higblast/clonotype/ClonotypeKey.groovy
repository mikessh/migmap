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

package com.antigenomics.higblast.clonotype

import com.antigenomics.higblast.genomic.Segment
import com.antigenomics.higblast.mapping.ReadMapping
import com.antigenomics.higblast.mutation.Mutation

class ClonotypeKey {
    final String cdr3nt
    final Segment vSegment, dSegment, jSegment
    final List<Mutation> mutations

    ClonotypeKey(ReadMapping readMapping) {
        this.cdr3nt = readMapping.cdr3nt
        this.vSegment = readMapping.mapping.vSegments[0]
        this.dSegment = readMapping.mapping.dSegments[0]
        this.jSegment = readMapping.mapping.jSegments[0]
        this.mutations = readMapping.mapping.mutations
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        ClonotypeKey that = (ClonotypeKey) o

        cdr3nt == that.cdr3nt &&
                jSegment == that.jSegment &&
                mutations == that.mutations &&
                vSegment == that.vSegment
    }

    @Override
    int hashCode() {
        int result
        result = cdr3nt.hashCode()
        result = 31 * result + vSegment.hashCode()
        result = 31 * result + jSegment.hashCode()
        31 * result + mutations.hashCode()
    }
}
