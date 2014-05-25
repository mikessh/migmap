package igblastwrp.io

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
class SeqData {
    final List<Integer> readIds = new LinkedList<>()
    final int[] qual
    final int seqId

    SeqData(int seqId, String qual) {
        this.seqId = seqId
        this.qual = qual ? new int[qual.length()] : null
    }

    void append(int readId, String newQual) {
        readIds.add(readId)
        if (this.qual) {
            def qualCharArray = newQual.toCharArray()
            for (int i = 0; i < qualCharArray.length; i++)
                qual[i] += qualCharArray[i] - 33
        }
    }

    String computeQual() {
        if (qual) {
            for (int i = 0; i < qual.length; i++)
                qual[i] /= readIds.size()
            def qualCharArray = new char[qual.length]
            for (int i = 0; i < qual.length; i++)
                qualCharArray[i] = (char) (qual[i] + 33)
            return new String(qualCharArray)
        }
        return null
    }
}
