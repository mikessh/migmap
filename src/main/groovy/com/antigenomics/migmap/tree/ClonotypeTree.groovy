/*
 * Copyright 2014-2015 Mikhail Shugay
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.migmap.tree

import com.antigenomics.migmap.clonotype.Clonotype
import com.milaboratory.core.sequence.NucleotideAlphabet
import com.milaboratory.core.sequence.NucleotideSequence
import com.milaboratory.core.tree.SequenceTreeMap
import com.milaboratory.util.Factory

class ClonotypeTree {
    final SequenceTreeMap<NucleotideSequence, List<Clonotype>> cTree = new SequenceTreeMap<>(NucleotideSequence.ALPHABET)

    ClonotypeTree(List<Clonotype> clonotypes) {
        clonotypes.each {
            def lst = cTree.createIfAbsent(new NucleotideSequence(it.cdr3nt), new Factory<List<Clonotype>>() {
                @Override
                List<Clonotype> create() {
                    new ArrayList<Clonotype>()
                }
            })
            lst << it
        }
    }

    void collapseGermlineMutationsInCdr3() {

    }
}
