package igblastwrp.blast
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
class JRefSearcher {
    private static final char GAP = '-'
    private final Map<String, Integer> jRefMap = new HashMap<>()

    public JRefSearcher(String species, String gene, String chain, File jRefFile) {
        def reader
        try {
            reader = new FileReader(jRefFile)
        } catch (FileNotFoundException e) {
            println "IgBlastWrapper bundle missing. $e"
            System.exit(-1)
        }
        def line
        //human	TRA	TRAJ37*01	27	TGGCTCTGGCAACACAGGCAAACTAATCTTTGGGCAAGGGACAACTTTACAAGTAAAACCAG
        while ((line = reader.readLine()) != null) {
            def splitLine = line.split("\t")
            if (splitLine[0] == species && splitLine[1] == "$gene$chain".toString()) {
                jRefMap.put(splitLine[2], Integer.parseInt(splitLine[3]))
            }
        }
    }

    public int getJRefPoint(String jSegment, int qstart, String qseq, int sstart, String sseq) {
        int jRef = jRefMap[jSegment]
        _getJRefPoint(jRef, qstart, qseq, sstart, sseq)
    }

    private static int _getJRefPoint(int jRef, int qstart, String qseq, int sstart, String sseq) {
        if (sstart > jRef || jRef + 4 >= sseq.replaceAll("-", "").length() + sstart)
            return -1

        int jRefRel = jRef - sstart, jRefDelta = 0
        for (int i = 0; i < jRefRel; i++) {
            if (qseq[i] == GAP)
                jRefDelta--
            else if (sseq[i] == GAP) {
                jRefDelta--
                jRefRel++
            }
        }

        jRefRel + jRefDelta + qstart
    }
}
