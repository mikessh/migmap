package com.antigenomics.migmap.post.analysis

import com.antigenomics.migmap.clonotype.Clonotype
import com.antigenomics.migmap.mutation.Mutation
import groovy.transform.CompileStatic

@CompileStatic
class Edge {
    final Clonotype from, to
    final Collection<Mutation> difference

    Edge(Clonotype from, Clonotype to, Collection<Mutation> difference) {
        this.from = from
        this.to = to
        this.difference = difference
    }
}
