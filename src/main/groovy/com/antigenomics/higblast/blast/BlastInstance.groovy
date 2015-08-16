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

import cc.redberry.pipe.OutputPortCloseable
import cc.redberry.pipe.VoidProcessor
import com.antigenomics.higblast.io.Read
import com.antigenomics.higblast.mapping.ReadMapping

class BlastInstance implements VoidProcessor<Read>, OutputPortCloseable<ReadMapping> {
    final Process process
    final BufferedReader reader
    final PrintWriter writer
    final BlastParser parser

    private List<String> header = []

    protected BlastInstance(Process process, BlastParser parser) {
        this.process = process
        this.reader = process.in.newReader()
        this.writer = process.out.newPrintWriter()
        this.parser = parser
    }

    private Read getCurrentRead() {
        new Read(header[0], header[1], header[2])
    }

    private ReadMapping parse(String chunk) {
        new ReadMapping(parser.parse(chunk), currentRead)
    }

    @Override
    synchronized ReadMapping take() {
        if (header == null) {
            // finished
            return null
        }

        ReadMapping readMapping
        String chunk = "", line

        while (true) {
            line = reader.readLine()

            if (!line) {
                if (chunk.length() == 0) {
                    return null // empty input
                } else {
                    readMapping = parse(chunk) // parse last chunk
                }
                header = null // warn of EOF
                break
            } else if (line.startsWith("# IGBLASTN") && chunk.length() > 0) {
                // went to next chunk, extract header and parse
                // header of next chunk will be with next readLine()
                header = (chunk =~ /# Query: (.+)\\|(.+)\\|(.+)/)[0][1..3]
                readMapping = parse(chunk)
                break
            } else {
                chunk += line + "\n"
            }
        }

        readMapping
    }

    @Override
    void process(Read input) {
        if (input) {
            writer.println(input.header + "|" + input.seq + "|" + input.qual)
        } else {
            // no more reads
            writer.close()
        }
    }

    @Override
    void close() {
        reader.close()
        process.waitFor()

        if (process.exitValue()) {
            println "[ERROR] ${process.getErrorStream()}"
        }
    }
}
