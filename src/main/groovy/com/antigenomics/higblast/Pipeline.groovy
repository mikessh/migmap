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

package com.antigenomics.higblast

import com.antigenomics.higblast.blast.BlastInstance
import com.antigenomics.higblast.blast.BlastInstanceFactory
import com.antigenomics.higblast.io.Read
import com.antigenomics.higblast.mapping.ReadMapping

class Pipeline {
    final OutputPort<Read> input
    final BlastInstanceFactory blastInstanceFactory
    final InputPort<ReadMapping> output

    Pipeline(OutputPort<Read> input, BlastInstanceFactory blastInstanceFactory,
             InputPort<ReadMapping> output) {
        this.input = input
        this.blastInstanceFactory = blastInstanceFactory
        this.output = output
    }

    void run(int nThreads = RuntimeInfo.N_THREADS) {
        def threads = new Thread[2 * nThreads]
        def instances = new BlastInstance[nThreads]

        // read >> blast instance threads
        (0..<nThreads).each {
            def instance = blastInstanceFactory.create()
            instances[it] = instance
            threads[it] = new Thread(new Runnable() {
                @Override
                void run() {
                    def read
                    while ((read = input.take()) != null) {
                        instance.put(read)
                    }
                    instance.put(null) // finished
                }
            })
        }

        // NOTE: here blast instance acts like a buffer

        // blast instance >> output threads
        (0..<nThreads).each {
            def instance = instances[it]
            threads[nThreads + it] = new Thread(new Runnable() {
                @Override
                void run() {
                    def result
                    while ((result = instance.take()) != null) {
                        output.put(result)
                    }
                }
            })
        }

        threads.each { it.start() }
        threads.each { it.join() }

        output.close()
    }
}
