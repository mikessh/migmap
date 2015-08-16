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

package com.antigenomics.higblast.blast

import cc.redberry.pipe.VoidProcessor
import cc.redberry.pipe.VoidProcessorFactory
import com.antigenomics.higblast.RuntimeInfo
import com.antigenomics.higblast.genomic.SegmentDatabase
import com.antigenomics.higblast.io.Read

class BlastInstanceFactory implements VoidProcessorFactory<Read> {
    final List cmdLine
    final File dir
    final List<String> env
    final BlastParser parser
    final SegmentDatabase segmentDatabase

    BlastInstanceFactory(String dataBundlePath,
                         String species, Set<String> genes,
                         boolean allAlleles, boolean useKabat) {
        this.segmentDatabase = new SegmentDatabase(dataBundlePath, species, genes, allAlleles, useKabat)
        this.parser = new BlastParser(segmentDatabase)

        def OPTS = ["SEGM_OPT"  :
                            ["v", "d", "j"].collect { segment ->
                                "-germline_db_${segment.toUpperCase()} $segmentDatabase.dbPath/$segment"
                            },

                    //"AUX_OPT"   :
                    //        "-auxiliary_data $dataBundlePath/optional_file/${species}_gl.aux",

                    "ORG_OPT"   :
                            "-organism $species",

                    "DOMAIN_OPT":
                            "-domain_system $useKabat"
        ]

        def REPORT_OPT = ["-num_alignments_V 1", "-num_alignments_D 1", "-num_alignments_J 1"],
            OUTFMT_OPT = ["-outfmt", "7 qseqid qstart qseq sstart sseq"],
            OUTPUT_OPT = "-out -"

        this.cmdLine = [[RuntimeInfo.igBlast, OPTS.values(), REPORT_OPT, OUTPUT_OPT].
                                flatten().join(" ").split(" "),
                        OUTFMT_OPT
        ].flatten()

        this.env = ["IGDATA=\"$dataBundlePath\""]
        this.dir = new File(dataBundlePath)

        segmentDatabase.makeBlastDb()
    }

    @Override
    VoidProcessor<Read> create() {
        new BlastInstance(cmdLine.execute(env, dir), parser)
    }
}
