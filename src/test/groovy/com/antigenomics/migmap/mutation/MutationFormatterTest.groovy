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

import com.antigenomics.migmap.blast.Alignment
import org.junit.Test

class MutationFormatterTest {
    @Test
    void reverseTest() {
        def alignment = new Alignment(0,
                //0000000001  11111111122222222223333
                //1234567890  12345678901234567890123
                "AGATCGATCGA--CTGCTACGACTGCATGACTCAAT", 0,
                //0000000001111111111222222   2222333
                //1234567890123456789012345   6789012
                "AAATCGATCGAAACTGCTACGACTGC---ACTCAGT")

        def mutations = MutationExtractor.extract(alignment).mutations

        def sseq = alignment.sseq.replaceAll("-", ""), qseq = alignment.qseq.replaceAll("-", "")

        assert sseq == MutationFormatter.mutateBack(qseq, mutations)
    }
}
