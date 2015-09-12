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

import com.antigenomics.higblast.ExecutionUtil
import com.antigenomics.higblast.Util
import com.antigenomics.higblast.blast.BlastInstanceFactory
import com.antigenomics.higblast.io.Read
import com.antigenomics.higblast.mapping.RegionMarkup

class SegmentDatabase {
    private static final List<SegmentDatabase> DB_CACHE = new LinkedList<>()
    final String databaseTempPath
    final Map<String, Segment> segments = new HashMap<>()
    final boolean hasD
    final int vSegments, dSegments, jSegments, vSegmentsNoMarkup
    private int annotatedV =0

    static final SPECIES_ALIAS =
            ["human"        : "HomoSapiens",
             "mouse"        : "MusMusculus",
             "rat"          : "RattusNorvegicus",
             "rabbit"       : "OryctolagusCuniculus",
             "rhesus_monkey": "MacacaMulatta"]

    SegmentDatabase(String dataBundlePath,
                    String species, List<String> genes,
                    boolean allAlleles = true) {
        genes = genes.unique()

        def speciesAlias = SPECIES_ALIAS[species]

        boolean hasD = false

        int vSegments = 0, dSegments = 0, jSegments = 0, vSegmentsNoMarkup = 0

        Util.getStream("segments.txt", true).splitEachLine("[\t ]+") { splitLine ->
            if (!splitLine[0].startsWith("#") &&
                    splitLine[0].startsWith(speciesAlias) &&
                    genes.contains(splitLine[1])) {

                def segmentName = splitLine[3]

                if (allAlleles || segmentName.endsWith("*01")) {
                    def seq = splitLine[5],
                        referencePoint = splitLine[4].toInteger(), segmentTypeStr = splitLine[2]

                    assert !segments.containsKey(segmentName)

                    if (segmentTypeStr.startsWith("V")) {
                        segments.put(segmentName, new Segment(SegmentType.V, segmentName, seq, referencePoint))
                        vSegments++
                    } else if (segmentTypeStr.startsWith("D")) {
                        segments.put(segmentName, new Segment(SegmentType.D, segmentName, seq, referencePoint))
                        hasD = true
                        dSegments++
                    } else if (segmentTypeStr.startsWith("J")) {
                        segments.put(segmentName, new Segment(SegmentType.J, segmentName, seq, referencePoint))
                        jSegments++
                    }
                }
            }
        }

        Util.report("Loaded database for $speciesAlias ${genes.join(",")} gene(s): " +
                "$vSegments Variable, $dSegments Diversity and $jSegments Joining segments.", 2)

        if (!hasD) {
            segments.put(Segment.DUMMY_D.name, Segment.DUMMY_D)
        }

        this.hasD = hasD
        this.databaseTempPath = dataBundlePath + "/database-" + UUID.randomUUID().toString()
        this.vSegments = vSegments
        this.dSegments = dSegments
        this.jSegments = jSegments
        this.vSegmentsNoMarkup = vSegmentsNoMarkup

        DB_CACHE.add(this)
    }

    void makeBlastDb() {
        Util.report("Creating temporary blast database at $databaseTempPath.", 2)

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
            "$ExecutionUtil.makeDb -parse_seqids -dbtype nucl -in $databaseTempPath/${it}.fa -out $databaseTempPath/$it".execute().waitFor()
        }
    }

    void annotateV(BlastInstanceFactory blastInstanceFactory) {
        Util.report("Annotating variable segments.", 2)

        def instance = blastInstanceFactory.create()

        def readThread = new Thread(new Runnable() {
            @Override
            void run() {
                def readMapping
                while ((readMapping = instance.take()) != null) {
                    def markup = readMapping.mapped ?
                            readMapping.mapping.regionMarkup : RegionMarkup.DUMMY
                    segments[readMapping.read.header].regionMarkup = markup
                    if (markup.complete) {
                        annotatedV++
                    }
                }
            }
        })

        readThread.start()

        segments.values()
                .findAll { it.type == SegmentType.V }
                .each {
            instance.put(new Read(it.name, it.sequence))
        }

        instance.put(null)

        readThread.join()

        Util.report("Fully annotated $annotatedV of $vSegments segments.", 2)
    }

    int getAnnotatedV() {
        return annotatedV
    }

    static void clearTemporaryFiles() {
        DB_CACHE.each {
            new File(it.databaseTempPath).deleteDir()
        }
    }
}
