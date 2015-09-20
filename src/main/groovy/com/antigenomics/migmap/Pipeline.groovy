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

package com.antigenomics.migmap

import com.antigenomics.migmap.blast.BlastInstance
import com.antigenomics.migmap.blast.BlastInstanceFactory
import com.antigenomics.migmap.io.InputPort
import com.antigenomics.migmap.io.OutputPort
import com.antigenomics.migmap.io.Read
import com.antigenomics.migmap.mapping.ReadMapping
import com.antigenomics.migmap.mapping.ReadMappingFilter

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
