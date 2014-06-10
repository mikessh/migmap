package igblastwrp.blast

import igblastwrp.Util
import igblastwrp.shm.Hypermutation

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
class ClonotypeData {
    final int level
    int nReads = 0, nEvents = 0
    final int[] qual1, qual2, qual3
    final Map<Hypermutation, Integer> hypermMap
    final Map<String, Integer> dSegmentVoter, vSegmentVoter, jSegmentVoter

    public ClonotypeData(String qual1, String qual2, String qual3, List<Hypermutation> hypermutations,
                         String vSegment, String dSegment, String jSegment, int nReads, int nEvents, int level) {
        this.level = level

        this.qual1 = (qual1 != Util.MY_NA) ? new byte[qual1.length()] : null
        this.qual2 = (qual2 != Util.MY_NA) ? new byte[qual2.length()] : null
        this.qual3 = (qual3 != Util.MY_NA) ? new byte[qual3.length()] : null

        this.hypermMap = new HashMap<>()
        this.vSegmentVoter = new HashMap<>()
        this.dSegmentVoter = new HashMap<>()
        this.jSegmentVoter = new HashMap<>()

        append(qual1, qual2, qual3, hypermutations, vSegment, dSegment, jSegment, nReads, nEvents)
    }

    void append(String qual1, String qual2, String qual3, List<Hypermutation> hypermutations,
                String vSegment, String dSegment, String jSegment, int nReads, int nEvents) {
        this.nEvents += nEvents
        this.nReads += nReads

        hypermutations.each {
            hypermMap.put(it, (hypermMap[it] ?: 0) + 1)
        }

        vSegmentVoter.put(vSegment, (vSegmentVoter[vSegment] ?: 0) + 1)
        dSegmentVoter.put(dSegment, (dSegmentVoter[dSegment] ?: 0) + 1)
        jSegmentVoter.put(jSegment, (jSegmentVoter[jSegment] ?: 0) + 1)

        if (qual1)
            if (qual1 != Util.MY_NA)
                for (int i = 0; i < qual1.length(); i++)
                    this.qual1[i] = this.qual1[i] + (int) (qual1.charAt(i) - 33) * nReads
        if (qual2)
            if (qual2 != Util.MY_NA)
                for (int i = 0; i < qual2.length(); i++)
                    this.qual2[i] = this.qual2[i] + (int) (qual2.charAt(i) - 33) * nReads
        if (qual3)
            if (qual3 != Util.MY_NA)
                for (int i = 0; i < qual3.length(); i++)
                    this.qual3[i] = this.qual3[i] + (int) (qual3.charAt(i) - 33) * nReads
    }

    byte summarizeQuality() {
        byte minQual = Byte.MAX_VALUE

        if (qual1)
            for (int i = 0; i < qual1.length; i++) {
                byte q = qual1[i] / nReads
                qual1[i] = q
                minQual = minQual > q ? q : minQual
            }
        if (qual2)
            for (int i = 0; i < qual2.length; i++) {
                byte q = qual2[i] / nReads
                qual2[i] = q
                minQual = minQual > q ? q : minQual
            }
        if (qual3)
            for (int i = 0; i < qual3.length; i++) {
                byte q = qual3[i] / nReads
                qual3[i] = q
                minQual = minQual > q ? q : minQual
            }

        minQual
    }

    String toString() {
        StringBuilder sb = new StringBuilder()

        sb.append(vSegmentVoter.max { it.value }.key).append('\t')
        sb.append(dSegmentVoter.max { it.value }.key).append('\t')
        sb.append(jSegmentVoter.max { it.value }.key).append('\t')

        if (qual1 && level > 0)
            qual1.each { sb.append((char) (it + 33)) }
        else
            sb.append(Util.MY_NA)
        sb.append('\t')

        if (qual2 && level > 0)
            qual2.each { sb.append((char) (it + 33)) }
        else
            sb.append(Util.MY_NA)
        sb.append('\t')

        if (qual3)
            qual3.each { sb.append((char) (it + 33)) }
        else
            sb.append(Util.MY_NA)
        sb.append('\t')

        if (hypermMap.size() == 0)
            sb.append(Util.MY_NA)
        else
            sb.append(hypermMap.sort { it.key.pos }.collect { it.value + "," + it.key.toString() }.join("|"))

        sb.toString()
    }

    final static String VALUE_HEADER = "vSegment\tdSegment\tjSegment\tcdr1q\tcdr2q\tcdr3q\tmutations"
}
