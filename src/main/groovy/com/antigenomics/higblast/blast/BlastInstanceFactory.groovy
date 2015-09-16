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

package com.antigenomics.higblast.blast

import com.antigenomics.higblast.ExecutionUtil
import com.antigenomics.higblast.genomic.SegmentDatabase

class BlastInstanceFactory {
    final List cmdLine
    final File dir
    final List<String> env
    final BlastParser parser
    final SegmentDatabase segmentDatabase

    BlastInstanceFactory(String dataBundlePath,
                         String species, List<String> genes,
                         boolean allAlleles, boolean useKabat) {
        this.segmentDatabase = new SegmentDatabase(dataBundlePath, species, genes, allAlleles)
        this.parser = new BlastParser(segmentDatabase)

        def seqtype = genes[0].startsWith("TR") ? "TCR" : "Ig"

        def OPTS = ["SEGM_OPT"  :
                            ["v", "d", "j"].collect { segment ->
                                "-germline_db_${segment.toUpperCase()} ${new File(segmentDatabase.databaseTempPath).absolutePath}/$segment"
                            },

                    "AUX_OPT"   :
                            "-auxiliary_data $dataBundlePath/optional_file/${species}_gl.aux",

                    "ORG_OPT"   :
                            "-organism $species -ig_seqtype $seqtype",

                    "DOMAIN_OPT":
                            "-domain_system ${useKabat ? "kabat" : "imgt"}"
        ]

        def REPORT_OPT = ["-num_alignments_V 3", "-num_alignments_D 3", "-num_alignments_J 3"],
            OUTFMT_OPT = ["-outfmt", "7 sseqid qstart qseq sstart sseq"],
            OUTPUT_OPT = "-out -"

        this.cmdLine = [[ExecutionUtil.igBlast, OPTS.values(), REPORT_OPT, OUTPUT_OPT].
                                flatten().join(" ").split(" "),
                        OUTFMT_OPT
        ].flatten()

        this.env = ["IGDATA=\"$dataBundlePath\""]
        this.dir = new File(dataBundlePath) // <- this is most crucial here. Blast will not find internal_data otherwise

        segmentDatabase.makeBlastDb()
    }

    void annotateV() {
        segmentDatabase.annotateV(this)
    }

    BlastInstance create() {
        new BlastInstance(cmdLine.execute(env, dir), parser)
    }
}
