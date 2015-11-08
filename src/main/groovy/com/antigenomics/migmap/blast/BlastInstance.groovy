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

package com.antigenomics.migmap.blast

import com.antigenomics.migmap.Util
import com.antigenomics.migmap.genomic.SegmentDatabase
import com.antigenomics.migmap.io.InputPort
import com.antigenomics.migmap.io.OutputPort
import com.antigenomics.migmap.io.Read
import com.antigenomics.migmap.mapping.ReadMapping
import groovy.transform.CompileStatic

import java.util.regex.Pattern

class BlastInstance implements OutputPort<ReadMapping>, InputPort<Read> {
    final Queue<Read> reads = new LinkedList<>()
    final Process proc
    final BufferedReader reader
    final PrintWriter writer
    final BlastParser parser
    final SegmentDatabase segmentDatabase

    protected boolean last = false

    protected BlastInstance(Process proc, BlastParser parser, SegmentDatabase segmentDatabase) {
        this.proc = proc
        this.reader = proc.in.newReader()
        this.writer = proc.out.newPrintWriter()
        this.parser = parser
        this.segmentDatabase = segmentDatabase
    }

    static String getHeader(String chunk) {
        (chunk =~ /# Query: (.+)/)[0][1]
    }

    @CompileStatic
    String nextChunk() {
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
    @CompileStatic
    ReadMapping take() {
        def chunk = nextChunk(), read

        if (chunk && (read = reads.poll())) { // second condition protects from empty input
            return new ReadMapping(parser.parse(chunk), read)
        } else {
            close()
            return null
        }
    }

    static final Pattern ALLOWED_BASES = ~/^[ATGCUatgcuRrYySsWwKkMmBbDdHhVvNn]+$/,
                         AMBIGOUS_BASE = ~/[^ATGCUatgcu]/;

    @CompileStatic
    static boolean isGoodBlastSeq(String seq) {
        if (seq.length() == 0 || !seq =~ ALLOWED_BASES) {
            return false
        }

        int goodBases = seq.length() - seq.replaceAll(AMBIGOUS_BASE, "").length()

        if (goodBases / (double) seq.length() >= 0.7) {
            return false
        }

        true
    }

    @Override
    @CompileStatic
    void put(Read input) {
        if (input) { // should contain at least one non-ambiguous base, otherwise BLAST crashes
            if (isGoodBlastSeq(input.seq)) {
                reads.add(input)
                writer.println(">" + input.header + "\n" + input.seq)
            }
        } else {
            // no more reads
            writer.close()
        }
    }

    @Override
    void close() {
        proc.waitFor()

        if (proc.exitValue() > 0) {
            Util.error("IgBlast failed with CODE${proc.exitValue()}:\n${proc.getErrorStream()}", 2)
        }
    }
}
