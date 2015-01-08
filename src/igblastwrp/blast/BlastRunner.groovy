package igblastwrp.blast

import igblastwrp.shm.SHMExtractor

import java.util.concurrent.ConcurrentHashMap

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
class BlastRunner implements Runnable {
    final static boolean WIN = System.properties['os.name'].toLowerCase().contains('windows')
    final List cmdLine
    final List env
    final File dir
    final BlastProcessor processor
    final ConcurrentHashMap<String, Clonotype> clonotypeMap

    public BlastRunner(int THREADS, String path, String species, String gene, String chain,
                       boolean allAlleles, String inputFileName,
                       ConcurrentHashMap<String, Clonotype> clonotypeMap) {
        def IGBLAST_EXECUTABLE = "$path/bin/igblastn",
            IGBLAST_DATA = "$path/data",
            IGBLAST_DB_PATH = "$IGBLAST_DATA/database",
            IGBLAST_OPT_PATH = "$IGBLAST_DATA/optional_file"

        if (WIN)
            IGBLAST_EXECUTABLE += ".exe"

        def seqtype = gene == "TR" ? "TCR" : "Ig"

        def OPTS = ["SEGM_OPT"  :
                            ["V", "J", "D"].collect { segment ->
                                "-germline_db_${segment} $IGBLAST_DB_PATH/${species}_${gene}_${chain}${allAlleles ? "_all" : ""}_${segment}"
                            },

                    "AUX_OPT"   :
                            "-auxiliary_data $IGBLAST_OPT_PATH/${species}_gl.aux",

                    "ORG_OPT"   :
                            "-organism $species -ig_seqtype $seqtype",

                    "THREAD_OPT":
                            "-num_threads $THREADS",

                    "DOMAIN_OPT":
                            "-domain_system imgt"
        ]

        def REPORT_OPT = ["-num_alignments_V 1", "-num_alignments_D 1", "-num_alignments_J 1"],
            OUTFMT_TRICK = ["-outfmt", "7 qseqid qstart qseq sstart sseq"],
            OUTPUT_OPT = "-out -",
            INPUT_OPT = "-query $inputFileName"

        cmdLine = [[IGBLAST_EXECUTABLE, OPTS.values(), REPORT_OPT, OUTPUT_OPT, INPUT_OPT].
                           flatten().join(" ").split(" "),
                   OUTFMT_TRICK
        ].flatten()

        env = ["IGDATA=\"$IGBLAST_DATA\""]
        dir = new File(IGBLAST_DATA)

        def jRefSearcher = new JRefSearcher(species, gene, chain, new File("$IGBLAST_DATA/jref.txt"))
        def shmExtractor = new SHMExtractor(
                "$IGBLAST_DB_PATH/${species}_${gene}_${chain}${allAlleles ? "_all" : ""}_V.fa",
                "$IGBLAST_DATA/internal_data/$species/${species}.ndm.imgt"
        )
        this.processor = new BlastProcessor(chain, jRefSearcher, shmExtractor)
        this.clonotypeMap = clonotypeMap
    }

    public void runIdle() {
        //println "[${new Date()}] Executing ${cmdLine.join(" ")}"

        def proc = cmdLine.execute(env, dir)

        proc.in.eachLine { line ->
            println line
        }

        proc.out.close()
        proc.waitFor()

        if (proc.exitValue()) {
            println "[ERROR] ${proc.getErrorStream()}"
        }
    }

    public void run() {
        //println "[${new Date()}] Executing ${cmdLine.join(" ")}"

        def proc = cmdLine.execute(env, dir)
        def clonotype
        String queryName
        String chunk = ""
        proc.in.eachLine { line ->
            if (line.startsWith("# IGBLASTN") && chunk.length() > 0) {
                // # Query: MIG UMI:TAGTCTGTAGCA:161
                queryName = (chunk =~ /# Query: (.+)/)[0][1]

                clonotype = processor.processChunk(chunk)

                if (clonotype)
                    clonotypeMap.put(queryName, clonotype)

                chunk = ""
            } else
                chunk += "$line\n"
        }

        // last chunk
        clonotype = processor.processChunk(chunk)

        if (clonotype) {
            queryName = (chunk =~ /# Query: (.+)/)[0][1]
            clonotypeMap.put(queryName, clonotype)
        }

        proc.out.close()
        proc.waitFor()

        if (proc.exitValue()) {
            println "[ERROR] ${proc.getErrorStream()}"
        }
    }
}
