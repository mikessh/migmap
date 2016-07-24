package com.antigenomics.migmap.post.analysis

import com.antigenomics.migmap.clonotype.Clonotype
import groovy.transform.CompileStatic

@CompileStatic
class Counter {
    int numberOfClonotypes, count
    double frequency

    void update(Clonotype clonotype) {
        numberOfClonotypes++
        frequency += clonotype.freq
        count += clonotype.count
    }
}
