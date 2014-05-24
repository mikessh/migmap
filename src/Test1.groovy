import igblastwrp.BlastProcessor
import igblastwrp.Clonotype

/**
 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */



def proc = new BlastProcessor("human", "TR", "A")

def seq1 = "AATGCGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTATCATCGCCAGTGGTATCAACG" +
        "CAGAGTTAGTTCTGTTAGCATCTTTGGGACTGGGACCAGATTACAAGTCTTTCCAAATATCCAGAACCCTGACCCAAGACC"

def chunk1 = "# IGBLASTN 2.2.29+\n" +
        "# Query: MIG UMI:TAGTCTGTAGCA:161\n" +
        "# Database: /Users/mikesh/Programming/igblastwrp/igblast/database/human_TR_A_V /Users/mikesh/Programming/igblastwrp/igblast/database/human_TR_A_D /Users/mikesh/Programming/igblastwrp/igblast/database/human_TR_A_J\n" +
        "# Domain classification requested: imgt\n" +
        "\n" +
        "# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.\n" +
        "TRAV12-3*01\tN/A\tN/A\tVB\tNo\tN/A\tN/A\t+\n" +
        "\n" +
        "# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself\n" +
        "GCAGA\tN/A\tN/A\tN/A\tN/A\t\n" +
        "\n" +
        "# Hit table (the first field indicates the chain type of the hit)\n" +
        "# Fields: query id, q. start, query seq, s. start, subject seq\n" +
        "# 1 hits found\n" +
        "V\tMIG\t63\tCCAGTGGTATCAACGCAGA\t155\tCCAGTGGTAACAAAGAAGA\n"

proc.processChunk(chunk1)

def seq2 = "TCCTTCAGTCTCAAGATCTCAGACTCACAGCTGGGGGACACTGCGATGTATTTCTGTGCTTTCATGAAGCCCCGAGA" +
        "GGGGAACACCGACAAGCTCATCTTTGGGACTGGGACCAGATTACAAGTCTTTCCAAATATCCAGAACCCTGACCCAAGACA"

def chunk2 = "# IGBLASTN 2.2.29+\n" +
        "# Query: MIG UMI:GATTTAACTGGG:1466\n" +
        "# Database: /Users/mikesh/Programming/igblastwrp/igblast/database/human_TR_A_V /Users/mikesh/Programming/igblastwrp/igblast/database/human_TR_A_D /Users/mikesh/Programming/igblastwrp/igblast/database/human_TR_A_J\n" +
        "# Domain classification requested: imgt\n" +
        "\n" +
        "# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.\n" +
        "TRAV38-1*01\tTRAJ34*01\tVA\tNo\tIn-frame\tYes\t+\n" +
        "\n" +
        "# V-(D)-J junction details based on top germline gene matches (V end, V-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself\n" +
        "GAAGC\tCCCGAGAGGGG\tAACAC\t\n" +
        "\n" +
        "# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)\n" +
        "FR3-IMGT\t1\t57\t57\t57\t0\t0\t100\n" +
        "CDR3-IMGT (germline)\t58\t70\t13\t13\t0\t0\t100\n" +
        "Total\tN/A\tN/A\t70\t70\t0\t0\t100\n" +
        "\n" +
        "# Hit table (the first field indicates the chain type of the hit)\n" +
        "# Fields: query id, q. start, query seq, s. start, subject seq\n" +
        "# 2 hits found\n" +
        "V\tMIG\t1\tTCCTTCAGTCTCAAGATCTCAGACTCACAGCTGGGGGACACTGCGATGTATTTCTGTGCTTTCATGAAGC\t220\tTCCTTCAGTCTCAAGATCTCAGACTCACAGCTGGGGGACACTGCGATGTATTTCTGTGCTTTCATGAAGC\n" +
        "J\tMIG\t82\tAACACCGACAAGCTCATCTTTGGGACTGGGACCAGATTACAAGTCTTTCCAA\t7\tAACACCGACAAGCTCATCTTTGGGACTGGGACCAGATTACAAGTCTTTCCAA\n"

Clonotype clonotype = proc.processChunk(chunk2)
println Clonotype.HEADER
println clonotype.generateEntry(seq2, null)
