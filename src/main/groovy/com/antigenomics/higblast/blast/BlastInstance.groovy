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

import com.antigenomics.higblast.InputPort
import com.antigenomics.higblast.OutputPort
import com.antigenomics.higblast.io.Read
import com.antigenomics.higblast.mapping.ReadMapping

class BlastInstance implements OutputPort<ReadMapping>, InputPort<Read> {
    final Process proc
    final BufferedReader reader
    final PrintWriter writer
    final BlastParser parser

    protected boolean last = false

    protected BlastInstance(Process proc, BlastParser parser) {
        this.proc = proc
        this.reader = proc.in.newReader()
        this.writer = proc.out.newPrintWriter()
        this.parser = parser
    }

    static Read getRead(String chunk) {
        def header = (chunk =~ /# Query: (.+)\|(.+)\|(.+)/)[0]
        new Read(header[1], header[2], header[3])
    }

    ReadMapping parse(String chunk) {
        new ReadMapping(parser.parse(chunk), getRead(chunk))
    }

    synchronized String nextChunk() {
        if (last) {
            // finished
            return null
        }

        String chunk = "", line

        while (true) {
            line = reader.readLine()

            if (line == null) {
                if (chunk.length() == 0) {
                    return null // empty input
                }
                last = true // warn of EOF
                break // went to EOF
            } else if (line.startsWith("# IGBLASTN") && chunk.length() > 0) {
                break // went to the beginning of next chunk
            } else {
                chunk += line + "\n"
            }
        }

        chunk
    }

    @Override
    ReadMapping take() {
        String chunk = nextChunk()

        chunk ? parse(chunk) : null
    }

    @Override
    void put(Read input) {
        if (input) {
            writer.println(">" + input.header + "|" + input.seq + "|" + input.qual)
            writer.println(input.seq)
        } else {
            // no more reads
            writer.close()
        }
    }

    @Override
    void close() {
        proc.waitFor()

        if (proc.exitValue() > 0) {
            println "[ERROR] code=${proc.exitValue()} ${proc.getErrorStream()}"
        }
    }
}
