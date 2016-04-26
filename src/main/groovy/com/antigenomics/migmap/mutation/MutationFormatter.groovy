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

package com.antigenomics.migmap.mutation

import com.antigenomics.migmap.Util
import groovy.transform.CompileStatic

@CompileStatic
class MutationFormatter {
    static
    final String OUTPUT_HEADER_NT = SubRegion.REGION_LIST.collect { SubRegion it -> "mutations.nt." + it }.join("\t"),
                 OUTPUT_HEADER_AA = SubRegion.REGION_LIST.collect { SubRegion it -> "mutations.aa." + it }.join("\t")

    static String toStringNT(List<Mutation> mutations) {
        // assume mutations are sorted

        def mutationStrings = [""] * SubRegion.REGION_LIST.length

        mutations.each {
            int order = it.subRegion.order

            if (mutationStrings[order].length() > 0)
                mutationStrings[order] += ","

            mutationStrings[order] += it.toString()
        }

        mutationStrings.join("\t")
    }

    static String toStringAA(List<Mutation> mutations, String rawQuery, int vStartInRef, int vStartInQuery) {
        String ref = Util.translateLinear('N' * vStartInRef + mutateBack(rawQuery, mutations).substring(vStartInQuery)),
               query = Util.translateLinear('N' * vStartInRef + rawQuery.substring(vStartInQuery))

        // assume mutations are sorted

        def mutationStrings = [""] * SubRegion.REGION_LIST.length

        mutations.each {
            int order = it.subRegion.order

            if (mutationStrings[order].length() > 0)
                mutationStrings[order] += ","

            int startInQuery = (int) ((it.startInRead - vStartInQuery + vStartInRef) / 3),
                endInQuery = (int) ((it.endInRead - vStartInQuery + vStartInRef - 1) / 3),
                startInRef = (int) (it.start / 3),
                endInRef = (int) ((it.end - 1) / 3)

            mutationStrings[order] += it.type.shortName + ((int) (it.pos / 3)).toString() + ":" +
                    (it.type == MutationType.Insertion ? "" : ref[startInRef..endInRef]) +
                    (it.type == MutationType.Substitution ? ">" : "") +
                    (it.type == MutationType.Deletion ? "" : query[startInQuery..endInQuery])
        }

        mutationStrings.join("\t")
    }

    static String mutateBack(String readSeq, List<Mutation> mutations) {
        List<String> mutatedSeq = readSeq.toCharArray().collect { it.toString() }

        mutations.each { mut ->
            switch (mut.type) {
                case MutationType.Substitution:
                    mutatedSeq[mut.startInRead] = mut.ntFrom
                    break
                case MutationType.Deletion:
                    mutatedSeq[mut.endInRead] = mut.ntFrom + mutatedSeq[mut.endInRead]
                    break
                case MutationType.Insertion:
                    (mut.startInRead..<mut.endInRead).each { int it -> mutatedSeq[it] = "" }
                    break
            }
        }

        mutatedSeq.join("")
    }
}
