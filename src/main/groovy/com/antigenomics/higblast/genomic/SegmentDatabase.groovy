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

import com.antigenomics.higblast.RuntimeInfo
import com.antigenomics.higblast.Util

class SegmentDatabase {
    final String databaseTempPath
    final Map<String, Segment> segments = new HashMap<>()
    final boolean hasD
    final int vSegments, dSegments, jSegments, vSegmentsNoMarkup

    static final SPECIES_ALIAS =
            ["human"        : "HomoSapiens",
             "mouse"        : "MusMusculus",
             "rat"          : "RattusNorvegicus",
             "rabbit"       : "OryctolagusCuniculus",
             "rhesus_monkey": "MacacaMulatta"]

    SegmentDatabase(String dataBundlePath,
                    String species, List<String> genes,
                    boolean allAlleles, boolean useKabat) {
        genes = genes.unique()

        def speciesAlias = SPECIES_ALIAS[species]

        def markupMap = new HashMap<String, int[]>()
        def markupPath = dataBundlePath + "/internal_data/$species/${species}.ndm.${useKabat ? "kabat" : "imgt"}"

        new File(markupPath).splitEachLine("[\t ]+") { splitLine ->
            markupMap.put(splitLine[0],
                    [splitLine[3].toInteger() - 1,
                     splitLine[4].toInteger(),
                     splitLine[7].toInteger() - 1,
                     splitLine[8].toInteger(),
                     splitLine[10].toInteger()] as int[]
            )
        }

        boolean hasD = false

        int vSegments = 0, dSegments = 0, jSegments = 0, vSegmentsNoMarkup = 0

        Util.getStream("segments.txt", true).splitEachLine("[\t ]+") { splitLine ->
            if (!splitLine[0].startsWith("#") &&
                    splitLine[0].startsWith(speciesAlias) &&
                    genes.contains(splitLine[1])) {

                def segmentName = splitLine[3]

                if (allAlleles || segmentName.endsWith("*01")) {
                    assert !segments.containsKey(segmentName)

                    if (splitLine[2].startsWith("V")) {
                        def markup = markupMap[segmentName]
                        if (markup) {
                            int referencePoint = splitLine[4].toInteger()
                            assert referencePoint == markup[4] // check if two reference files are concordant
                            segments.put(segmentName, new VSegment(segmentName, splitLine[5], referencePoint,
                                    markup[0], markup[1], markup[2], markup[3]))
                            vSegments++
                        } else {
                            vSegmentsNoMarkup++
                            Util.report("[WARNING] No markup for $segmentName V segment.", 3)
                        }
                    } else if (splitLine[2].startsWith("D")) {
                        segments.put(segmentName, new DSegment(segmentName, splitLine[5]))
                        hasD = true
                        dSegments++
                    } else if (splitLine[2].startsWith("J")) {
                        segments.put(segmentName, new JSegment(segmentName, splitLine[5],
                                splitLine[4].toInteger()))
                        jSegments++
                    }
                }
            }
        }

        Util.report("Loaded database. #v=$vSegments,#d=$dSegments,#j=$jSegments.", 2)

        if (!hasD) {
            segments.put(DSegment.DUMMY.name, DSegment.DUMMY)
        }

        this.hasD = hasD
        this.databaseTempPath = dataBundlePath + "/database-" + UUID.randomUUID().toString()
        this.vSegments = vSegments
        this.dSegments = dSegments
        this.jSegments = jSegments
        this.vSegmentsNoMarkup = vSegmentsNoMarkup

    }

    void makeBlastDb() {
        new File(databaseTempPath).mkdir()

        new File("$databaseTempPath/v.fa").withPrintWriter { pwV ->
            new File("$databaseTempPath/d.fa").withPrintWriter { pwD ->
                new File("$databaseTempPath/j.fa").withPrintWriter { pwJ ->
                    segments.values().each {
                        switch (it.type) {
                            case SegmentType.V:
                                pwV.println(it.toFastaString())
                                break

                            case SegmentType.J:
                                pwJ.println(it.toFastaString())
                                break

                            case SegmentType.D:
                                pwD.println(it.toFastaString())
                                break
                        }
                    }
                }
            }
        }

        ["v", "d", "j"].each {
            "$RuntimeInfo.makeDb -parse_seqids -dbtype nucl -in $databaseTempPath/${it}.fa -out $databaseTempPath/$it".execute().waitFor()
        }
    }

    void clearBlastDb() {
        new File(databaseTempPath).deleteDir()
    }
}
