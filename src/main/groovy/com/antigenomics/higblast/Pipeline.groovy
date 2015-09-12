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
import com.antigenomics.higblast.io.InputPort
import com.antigenomics.higblast.io.OutputPort
import com.antigenomics.higblast.io.Read
import com.antigenomics.higblast.mapping.ReadMapping
import com.antigenomics.higblast.mapping.ReadMappingFilter

import java.util.concurrent.atomic.AtomicLong

class Pipeline {
    private final AtomicLong inputCounter = new AtomicLong()
    final long limit
    final int nThreads
    final OutputPort<Read> input
    final BlastInstanceFactory blastInstanceFactory
    final InputPort<ReadMapping> output
    final ReadMappingFilter readMappingFilter

    Pipeline(OutputPort<Read> input, BlastInstanceFactory blastInstanceFactory,
             InputPort<ReadMapping> output, ReadMappingFilter readMappingFilter,
             long limit = -1, int nThreads = Util.N_THREADS) {
        this.input = input
        this.blastInstanceFactory = blastInstanceFactory
        this.output = output
        this.readMappingFilter = readMappingFilter
        this.limit = limit < 0 ? Long.MAX_VALUE : limit
        this.nThreads = nThreads
    }

    void run() {
        def threads = new Thread[2 * nThreads]
        def instances = new BlastInstance[nThreads]

        Util.report("Started analysis", 2)

        // read >> blast instance threads
        (0..<nThreads).each {
            def instance = blastInstanceFactory.create()
            instances[it] = instance
            threads[it] = new Thread(new Runnable() {
                @Override
                void run() {
                    def read
                    while (((read = input.take()) != null) &&
                            (inputCounter.incrementAndGet() <= limit)
                    ) {
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
                        if (readMappingFilter.pass(result)) {
                            output.put(result)
                        }
                    }
                }
            })
        }

        def reporter = new Thread(new Runnable() {
            @Override
            void run() {
                try {
                    while (!Thread.currentThread().isInterrupted()) {
                        Util.report("Loaded $inputCount reads. " +
                                (readMappingFilter.total > 0L ? readMappingFilter.toProgressString() : ""), 2)
                        Thread.sleep(10000)
                    }
                } catch (InterruptedException e) {

                }
            }
        })

        threads.each { it.start() }

        reporter.daemon = true
        reporter.start()

        threads.each { it.join() }

        reporter.interrupt()

        Util.report("Finished analysis. ${readMappingFilter.toProgressString()}", 2)

        // Close all ports

        readMappingFilter.unmappedInputPort.close()
        output.close()
    }

    final Thread.UncaughtExceptionHandler h = new Thread.UncaughtExceptionHandler() {
        void uncaughtException(Thread t, Throwable e) {
            throw new RuntimeException("Error in pipeline", e)
        }
    }

    long getInputCount() {
        inputCounter.get()
    }
}
