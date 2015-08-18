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

package com.antigenomics.higblast.mapping

import com.antigenomics.higblast.genomic.DSegment
import com.antigenomics.higblast.genomic.JSegment
import com.antigenomics.higblast.genomic.VSegment
import com.antigenomics.higblast.mutation.Mutation

class Mapping {
    final List<VSegment> vSegments
    final List<DSegment> dSegments
    final List<JSegment> jSegments
    final boolean rc, complete, hasCdr3, inFrame, noStop, hasD
    final Cdr3Markup cdr3Markup
    final RegionMarkup regionMarkup
    final List<Mutation> mutations

    Mapping(List<VSegment> vSegments, List<DSegment> dSegments, List<JSegment> jSegments,
            RegionMarkup regionMarkup, Cdr3Markup cdr3Markup,
            boolean rc, boolean complete, boolean hasCdr3, boolean inFrame, boolean noStop, hasD,
            List<Mutation> mutations) {
        // needed to extract CDR3 sequence
        this.regionMarkup = regionMarkup
        this.rc = rc

        // main data
        this.vSegments = vSegments
        this.dSegments = dSegments
        this.jSegments = jSegments

        this.cdr3Markup = cdr3Markup
        this.mutations = mutations

        // misc
        this.complete = complete
        this.hasCdr3 = hasCdr3
        this.inFrame = inFrame
        this.noStop = noStop
        this.hasD = hasD
    }
}
