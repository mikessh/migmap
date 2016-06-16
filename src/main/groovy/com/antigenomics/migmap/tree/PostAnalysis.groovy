package com.antigenomics.migmap.tree

import com.antigenomics.migmap.clonotype.Clonotype
import com.antigenomics.migmap.mutation.Mutation

class PostAnalysis {
    final List<Clonotype> sample

    PostAnalysis(List<Clonotype> clonotypes) {
        this.sample = clonotypes
    }

    void generateHypermutationTable(String fileName) {
        def mutationMap = new HashMap<Mutation, Counter>()

        sample.each { Clonotype clonotype ->
            clonotype.mutations.each { Mutation mutation ->
                Counter counter
                mutationMap.put(mutation, counter = (mutationMap[mutation] ?: new Counter()))
                counter.update(clonotype)
            }
        }

        new File(fileName).withPrintWriter { pw ->
            pw.println("segment\tsegment.name\tregion\tmutation.type\t" +
                    "pos.nt\tfrom.nt\tto.nt\tpos.aa\tfrom.aa\tto.aa\t" +
                    "count.clonotypes\tcount.reads\tcount.freq")

            mutationMap.each {
                def mutation = it.key, counter = it.value
                pw.println([mutation.parent.gene,
                            mutation.parent.name,
                            mutation.subRegion,
                            mutation.type,
                            mutation.pos,
                            mutation.ntFrom,
                            mutation.ntTo,
                            mutation.aaPos,
                            mutation.aaFrom,
                            mutation.aaTo,
                            counter.numberOfClonotypes,
                            counter.count,
                            counter.frequency
                ])
            }
        }
    }


}
