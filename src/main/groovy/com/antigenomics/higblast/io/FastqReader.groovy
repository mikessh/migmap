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

package com.antigenomics.higblast.io

import cc.redberry.pipe.OutputPortCloseable
import com.antigenomics.higblast.Util

class FastqReader implements OutputPortCloseable<Read> {
    final BufferedReader reader

    FastqReader(String fileName) {
        reader = Util.getReader(fileName)
    }

    @Override
    void close() {
        reader.close()
    }

    @Override
    synchronized Read take() {
        def header = reader.readLine()
        if (!header) {
            return null
        }

        def seq = reader.readLine()
        reader.readLine()
        def qual = reader.readLine()

        new Read(header, seq, qual)
    }
}
