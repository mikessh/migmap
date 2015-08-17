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
    final String dbPath
    final Map<String, Segment> segments = new HashMap<>()
    final boolean hasD

    static final SPECIES_ALIAS =
            ["human"        : "HomoSapiens",
             "mouse"        : "MusMusculus",
             "rat"          : "RattusNorvegicus",
             "rabbit"       : "OryctolagusCuniculus",
             "rhesus_monkey": "MacacaMulatta"]

    SegmentDatabase(String dataBundlePath,
                    String species, Set<String> genes,
                    boolean allAlleles, boolean useKabat) {

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
                            segments.put(segmentName, new VSegment(name: segmentName, sequence: splitLine[5],
                                    cdr1start: markup[0], cdr1end: markup[1], cdr2start: markup[2], cdr2end: markup[3],
                                    referencePoint: referencePoint))
                            vSegments++
                        } else {
                            vSegmentsNoMarkup
                        }
                    } else if (splitLine[2].startsWith("D")) {
                        segments.put(segmentName, new DSegment(name: segmentName, sequence: splitLine[5]))
                        hasD = true
                        dSegments++
                    } else if (splitLine[2].startsWith("J")) {
                        segments.put(segmentName, new JSegment(name: segmentName, sequence: splitLine[5],
                                referencePoint: splitLine[4].toInteger()))
                        jSegments++
                    }
                }
            }
        }

        println "Loaded database. #v=$vSegments,#d=$dSegments,#j=$jSegments. No markup for $vSegmentsNoMarkup v"

        if (!hasD) {
            segments.put(".", new DSegment(name: ".", sequence: "GGGGGGGGGGGGGGG"))
        }

        this.hasD = hasD
        this.dbPath = dataBundlePath + "/database-" + UUID.randomUUID().toString()
    }

    void makeBlastDb() {
        new File(dbPath).mkdir()

        new File("$dbPath/v.fa").withPrintWriter { pwV ->
            new File("$dbPath/d.fa").withPrintWriter { pwD ->
                new File("$dbPath/j.fa").withPrintWriter { pwJ ->
                    segments.values().each {
                        if (it instanceof VSegment)
                            pwV.println(it.toFastaString())
                        else if (it instanceof JSegment)
                            pwJ.println(it.toFastaString())
                        else if (it instanceof DSegment)
                            pwD.println(it.toFastaString())
                    }
                }
            }
        }

        ["v", "d", "j"].each {
            "$RuntimeInfo.makeDb -parse_seqids -dbtype nucl -in $dbPath/${it}.fa -out $dbPath/$it".execute().waitFor()
        }
    }

    void clearBlastDb() {
        new File(dbPath).deleteDir()
    }
}
