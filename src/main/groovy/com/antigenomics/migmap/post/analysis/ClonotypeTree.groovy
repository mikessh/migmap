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

package com.antigenomics.migmap.post.analysis

import com.antigenomics.migmap.pipeline.Util
import com.antigenomics.migmap.clonotype.Clonotype
import com.antigenomics.migmap.mutation.Mutation
import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.NucleotideSequence
import com.milaboratory.core.tree.SequenceTreeMap
import com.milaboratory.core.tree.TreeSearchParameters
import com.milaboratory.util.Factory

class ClonotypeTree {
    static final TreeSearchParameters germOnlyDiffSearchParams = new TreeSearchParameters(5, 0, 0, 5)

    // Grouping by CDR3
    static Collection<List<Clonotype>> groupClonotypesByCdr3(Iterable<Clonotype> sample) {
        def cTree = new SequenceTreeMap<NucleotideSequence, List<Clonotype>>(NucleotideSequence.ALPHABET)
        def cdr3RepresentativeMap = new HashMap<NucleotideSequence, Clonotype>()

        Util.report("- Pre-group by CDR3")
        sample.each { Clonotype clonotype ->
            def cdr3Nt = new NucleotideSequence(clonotype.cdr3nt)
            cdr3RepresentativeMap.putIfAbsent(cdr3Nt, clonotype)
            def lst = cTree.createIfAbsent(cdr3Nt, new Factory<List<Clonotype>>() {
                @Override
                List<Clonotype> create() {
                    def list = new ArrayList<Clonotype>()
                    list
                }
            })
            lst << clonotype
        }

        Util.report("- Finding CDR3 differences explained by SHMs in germline region")
        def cdr3MatchWithoutGermlineNet = new HashMap<NucleotideSequence, CNode>()
        cdr3RepresentativeMap.each {
            def niter = cTree.getNeighborhoodIterator(it.key, germOnlyDiffSearchParams)
            def cNode
            cdr3MatchWithoutGermlineNet.put(it.key, cNode = new CNode(cdr3MatchWithoutGermlineNet, it.key))

            def nextClonotypeList
            while ((nextClonotypeList = niter.next()) != null) {
                def otherRepresentativeClone = nextClonotypeList.first()
                def otherCdr3Nt = new NucleotideSequence(otherRepresentativeClone.cdr3nt)

                if (!cdr3MatchWithoutGermlineNet.containsKey(otherCdr3Nt) && !it.key.equals(otherCdr3Nt)) {
                    def alignment = niter.currentAlignment

                    if (cdr3MatchWithoutGermlineMutations(it.value, otherRepresentativeClone,
                            otherCdr3Nt, alignment)) {
                        cNode.children.add(otherCdr3Nt)
                    }
                }
            }
        }

        Util.report("- Aggregating")
        cdr3MatchWithoutGermlineNet.each {
            if (!it.value.parent) {
                it.value.assignParent(it.key)
            }
        }

        def groupedMap = new HashMap<NucleotideSequence, List<Clonotype>>()

        cdr3MatchWithoutGermlineNet.values().each {
            def clonotypeList
            groupedMap.put(it.parent, clonotypeList = (groupedMap[it.parent] ?: new ArrayList<Clonotype>()))
            clonotypeList.addAll(cTree.get(it.tag))
        }

        groupedMap.values()
    }

    private static class CNode {
        final HashMap<NucleotideSequence, CNode> net
        final List<NucleotideSequence> children = new ArrayList<>()
        final NucleotideSequence tag
        NucleotideSequence parent

        CNode(HashMap<NucleotideSequence, CNode> net, NucleotideSequence tag) {
            this.net = net
            this.tag = tag
        }

        void assignParent(NucleotideSequence parent) {
            this.parent = parent
            children.each {
                net[it].assignParent(parent)
            }
        }
    }

    static boolean inGermline(int position, Clonotype clonotype) {
        clonotype.cdr3Markup.vEnd >= position || clonotype.cdr3Markup.jStart <= position ||
                (clonotype.cdr3Markup.dEnd >= position && clonotype.cdr3Markup.dStart <= position)
    }

    static boolean cdr3MatchWithoutGermlineMutations(Clonotype clonotype1, Clonotype clonotype2,
                                                     NucleotideSequence cdr3nt2,
                                                     Alignment<NucleotideSequence> aln) {
        def mut12 = aln.absoluteMutations, mut21 = aln.invert(cdr3nt2).absoluteMutations

        // in germline in at least one of sequences

        for (int i = 0; i < mut12.size(); i++) {
            if (!inGermline(mut12.getPositionByIndex(i), clonotype1) &&
                    !inGermline(mut21.getPositionByIndex(i), clonotype2))
                return false
        }
        true
    }

    // Graph

    static List<Edge> getEdges(List<Clonotype> clonotypes) {
        // Try connect clonotypes
        def edgeList = new ArrayList<Edge>()
        def parentNodeMap = new HashMap<Clonotype, Set<Clonotype>>()
        for (int i = 0; i < clonotypes.size(); i++) {
            for (int j = i + 1; j < clonotypes.size(); j++) {
                def edge = connect(clonotypes[i], clonotypes[j])
                if (edge) {
                    edgeList.add(edge)
                    def parentList
                    parentNodeMap.put(edge.to,
                            parentList = (parentNodeMap[edge.to] ?: new HashSet<>()))
                    parentList.add(edge.from)
                }
            }
        }

        // Remove redundancy
        edgeList.findAll { edge ->
            //
            //      2 <-- 3
            //      ^     ^
            //       \   /
            //         1
            //
            // remove 1->2
            // any parent(2) that is not 1 has 1 as parent
            !parentNodeMap[edge.to].any { it != edge.from && hasParent(parentNodeMap[it], edge.from) }
        }
    }

    private static boolean hasParent(Set<Clonotype> clonotypes, Clonotype clonotype) {
        clonotypes && clonotypes.contains(clonotype)
    }

    static Set<String> getIntersection(List<Mutation> mutations1, List<Mutation> mutations2) {
        def set = new HashSet<String>(mutations1*.toString())
        new HashSet<>(mutations2*.toString().findAll { set.contains(it.toString()) })
    }

    static Edge connect(Clonotype clonotype1, Clonotype clonotype2) {
        def mutations1 = clonotype1.mutations,
            mutations2 = clonotype2.mutations

        def intersection = mutations1.size() > mutations2.size() ?
                getIntersection(mutations1, mutations2) : getIntersection(mutations2, mutations1)

        intersection.size() == mutations1.size() ?
                new Edge(clonotype1, clonotype2, mutations2.findAll {
                    !intersection.contains(it)
                }) :
                (intersection.size() == mutations2.size() ?
                        new Edge(clonotype2, clonotype1, mutations1.findAll {
                            !intersection.contains(it)
                        }) :
                        null)
    }
}
