package com.antigenomics.higblast.io

import com.antigenomics.higblast.Util

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
class FastaReader {
    final BufferedReader reader
    String header = ""

    FastaReader(String fileName) {
        reader = Util.getReader(fileName)
    }

    Read next() {
        if (header == "") // first sequence
            header = reader.readLine()

        if (!header) // EOF
            return null

        def seq = ""
        while (true) {
            def line = reader.readLine()
            if (line == null) {
                // EOF next read, return last sequence
                def _header = header
                header = null
                return new Read(_header, seq, null)
            } else if (line.startsWith(">")) {
                // reset header and return current read
                def _header = header
                header = line
                return new Read(_header, seq, null)
            } else // sequence line
                seq += line
        }
    }
}
