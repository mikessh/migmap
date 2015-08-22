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

import java.util.concurrent.atomic.AtomicLong

class Pipeline {
    private final AtomicLong inputCounter = new AtomicLong(),
                             processedCounter = new AtomicLong(),
                             mappedCounter = new AtomicLong(),
                             cdr3FoundCounter = new AtomicLong()
    final long limit
    final int nThreads
    final OutputPort<Read> input
    final BlastInstanceFactory blastInstanceFactory
    final InputPort<ReadMapping> output

    Pipeline(OutputPort<Read> input, BlastInstanceFactory blastInstanceFactory,
             InputPort<ReadMapping> output, long limit = -1, int nThreads = Util.N_THREADS) {
        this.input = input
        this.blastInstanceFactory = blastInstanceFactory
        this.output = output
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
                        if (result.mapped) {
                            mappedCounter.incrementAndGet()
                            if (result.mapping.hasCdr3) {
                                cdr3FoundCounter.incrementAndGet()
                            }
                        }
                        processedCounter.incrementAndGet()
                        output.put(result)
                    }
                }
            })
        }

        def reporter = new Thread(new Runnable() {
            @Override
            void run() {
                try {
                    while (!Thread.currentThread().isInterrupted()) {
                        Util.report("Loaded $inputCount reads, processed $processedCount. " +
                                (processedCount > 0L ?
                                        ("Mapping efficincy ${(int) (mappedRatio * 10000) / 100}% " +
                                                "(with CDR3 ${(int) (cdr3FoundRatio * 10000) / 100}%)") :
                                        ""), 2)
                        Thread.sleep(10000)
                    }
                } catch (InterruptedException e) {

                }
            }
        })

        threads.each { it.start() }

        reporter.start()

        threads.each { it.join() }

        reporter.interrupt()

        Util.report("Finished analysis. Reads mapped $mappedCount (with CDR3 $cdr3FoundCount) of $processedCount", 2)

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

    long getProcessedCount() {
        processedCounter.get()
    }

    long getMappedCount() {
        mappedCounter.get()
    }

    long getCdr3FoundCount() {
        cdr3FoundCounter.get()
    }

    double getMappedRatio() {
        (double) mappedCount / processedCount
    }

    double getCdr3FoundRatio() {
        (double) cdr3FoundCount / processedCount
    }
}
