package com.antigenomics.migmap.tree

import com.antigenomics.migmap.clonotype.Clonotype
import com.antigenomics.migmap.genomic.Segment
import com.antigenomics.migmap.mutation.Mutation

class PostAnalysis {
    final List<Clonotype> sample

    PostAnalysis(List<Clonotype> clonotypes) {
        this.sample = clonotypes
    }

    void generateClonotypeTree(String fileNamePrefix) {
        def clonotypeTree = new ClonotypeTree(sample)
    }

    void generateHypermutationTable(String fileName) {
        def segmentMap = new HashMap<Segment, Counter>()

        sample.each { Clonotype clonotype ->
            Counter counter
            [clonotype.vSegment, clonotype.dSegment, clonotype.jSegment].findAll { !Segment.isDummy(it) }.each {
                segmentMap.put(it, counter = (segmentMap[it] ?: new Counter()))
                counter.update(clonotype)
            }
        }

        new File(fileName).withPrintWriter { pw ->
            pw.println("clonotype.id\tclonotype.v\tclonotype.j\tsegment\tsegment.name\tregion\tmutation.type\t" +
                    "pos.nt\tfrom.nt\tto.nt\tpos.aa\tfrom.aa\tto.aa\t" +
                    "count.clonotypes\tcount.reads\tcount.freq\t" +
                    "total.clonotypes\ttotal.count\ttotal.freq")
            sample.eachWithIndex { clonotype, ind ->
                clonotype.mutations.each { mutation ->
                    def segmentCounter = segmentMap[mutation.parent]
                    pw.println([ind,
                                clonotype.vSegment.name,
                                clonotype.jSegment.name,
                                mutation.parent.type,
                                mutation.parent.name,
                                mutation.subRegion,
                                mutation.type,
                                mutation.pos,
                                mutation.ntFrom,
                                mutation.ntTo,
                                mutation.aaPos,
                                mutation.aaFrom,
                                mutation.aaTo,
                                1,
                                clonotype.count,
                                clonotype.freq,
                                segmentCounter.numberOfClonotypes,
                                segmentCounter.count,
                                segmentCounter.frequency
                    ].join("\t"))
                }
            }
        }
    }
}
