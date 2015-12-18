/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */

package com.antigenomics.migmap.genomic

import com.antigenomics.migmap.ExecutionUtil
import com.antigenomics.migmap.Util
import com.antigenomics.migmap.blast.BlastInstanceFactory
import com.antigenomics.migmap.io.Read
import com.antigenomics.migmap.mapping.RegionMarkup
import groovy.transform.CompileStatic

@CompileStatic
class SegmentDatabase {
    private static final List<SegmentDatabase> DB_CACHE = new LinkedList<>()
    final String databaseTempPath
    final Set<String> genes = new HashSet<>()
    final Map<String, Segment> segments = new HashMap<>()
    final boolean hasD
    final int vSegments, dSegments, jSegments
    private int annotatedV = 0

    static final SPECIES_ALIAS =
            ["human"        : "HomoSapiens",
             "mouse"        : "MusMusculus",
             "rat"          : "RattusNorvegicus",
             "rabbit"       : "OryctolagusCuniculus",
             "rhesus_monkey": "MacacaMulatta"]

    SegmentDatabase(String dataBundlePath,
                    String species, List<String> genes,
                    boolean allAlleles = true,
                    String segmentsFilePath = null) {
        this.genes.addAll(genes)

        String speciesAlias = SPECIES_ALIAS[species]

        boolean hasD = false

        int vSegments = 0, dSegments = 0, jSegments = 0

        segmentsFilePath ? new File(segmentsFilePath) : Util.getStream("segments.txt", true)
                .splitEachLine("[\t ]+") { List<String> splitLine ->
            if (!splitLine[0].startsWith("#") &&
                    splitLine[0].startsWith(speciesAlias) &&
                    this.genes.contains(splitLine[1])) {

                def gene = splitLine[1], segmentName = splitLine[3]

                if (allAlleles || segmentName.endsWith("*01")) {
                    def seq = splitLine[5],
                        referencePoint = splitLine[4].toInteger(), segmentTypeStr = splitLine[2]

                    assert !segments.containsKey(segmentName)

                    if (segmentTypeStr.startsWith("V")) {
                        segments.put(segmentName, new Segment(this, SegmentType.V, gene, segmentName, seq, referencePoint))
                        vSegments++
                    } else if (segmentTypeStr.startsWith("D")) {
                        segments.put(segmentName, new Segment(this, SegmentType.D, gene, segmentName, seq, referencePoint))
                        hasD = true
                        dSegments++
                    } else if (segmentTypeStr.startsWith("J")) {
                        segments.put(segmentName, new Segment(this, SegmentType.J, gene, segmentName, seq, referencePoint))
                        jSegments++
                    }
                }
            }
        }

        Util.report("Loaded database for $speciesAlias ${this.genes.join(",")} gene(s): " +
                "$vSegments Variable, $dSegments Diversity and $jSegments Joining segments.", 2)

        if (!hasD) {
            segments.put(Segment.DUMMY_D.name, Segment.DUMMY_D)
        }

        this.hasD = hasD
        this.databaseTempPath = dataBundlePath + "/database-" + UUID.randomUUID().toString()
        this.vSegments = vSegments
        this.dSegments = dSegments
        this.jSegments = jSegments

        DB_CACHE.add(this)
    }

    void makeBlastDb() {
        Util.report("Creating temporary BLAST database $databaseTempPath.", 2)

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
