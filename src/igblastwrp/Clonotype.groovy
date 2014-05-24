package igblastwrp

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
class Clonotype {
    final String vSegment, dSegment, jSegment
    final int cdr1start, cdr1end, cdr2start, cdr2end, cdr3start, cdr3end

    Clonotype(String vSegment, String dSegment, String jSegment,
              int cdr1start, int cdr1end, int cdr2start, int cdr2end, int cdr3start, int cdr3end) {
        this.vSegment = vSegment
        this.dSegment = dSegment
        this.jSegment = jSegment
        this.cdr1start = cdr1start
        this.cdr1end = cdr1end
        this.cdr2start = cdr2start
        this.cdr2end = cdr2end
        this.cdr3start = cdr3start
        this.cdr3end = cdr3end
    }

    String generateEntry(String seq, String qual) {
        def cdr1nt = cdr1start >= 0 ? seq.substring(cdr1start, cdr1end) : "N/A",
            cdr2nt = cdr2start >= 0 ? seq.substring(cdr2start, cdr2end) : "N/A",
            cdr3nt = cdr3start >= 0 ? seq.substring(cdr3start, cdr3end) : "N/A"

        def cdr1q = cdr1start >= 0 && qual ? qual.substring(cdr1start, cdr1end) : "N/A",
            cdr2q = cdr2start >= 0 && qual ? qual.substring(cdr2start, cdr2end) : "N/A",
            cdr3q = cdr3start >= 0 && qual ? qual.substring(cdr3start, cdr3end) : "N/A"

        def cdr1aa = cdr1start >= 0 ? translate(cdr1nt) : "N/A",
            cdr2aa = cdr2start >= 0 ? translate(cdr2nt) : "N/A",
            cdr3aa = cdr3start >= 0 ? translate(cdr3nt) : "N/A"

        [vSegment, dSegment, jSegment, cdr1nt, cdr2nt, cdr3nt, cdr1q, cdr2q, cdr3q, cdr1aa, cdr2aa, cdr3aa].join("\t")
    }


    final
    static HEADER = "v_segment\td_segment\tj_segment\tcdr1nt\tcdr2nt\tcdr3nt\tcdr1q\tcdr2q\tcdr3q\tcdr1aa\tcdr2aa\tcdr3aa",
           HEADER_RAW = "v_segment\td_segment\tj_segment\tcdr1start\tcdr1end\tcdr2start\tcdr2end\tcdr3start\tcdr3end"

    static String codon2aa(String codon) {
        String codonUpper = codon.toUpperCase()
        switch (codonUpper) {
            case 'TTT': return 'F'
            case 'TTC': return 'F'
            case 'TTA': return 'L'
            case 'TTG': return 'L'
            case 'TCT': return 'S'
            case 'TCC': return 'S'
            case 'TCA': return 'S'
            case 'TCG': return 'S'
            case 'TAT': return 'Y'
            case 'TAC': return 'Y'
            case 'TAA': return '*'
            case 'TAG': return '*'
            case 'TGT': return 'C'
            case 'TGC': return 'C'
            case 'TGA': return '*'
            case 'TGG': return 'W'
            case 'CTT': return 'L'
            case 'CTC': return 'L'
            case 'CTA': return 'L'
            case 'CTG': return 'L'
            case 'CCT': return 'P'
            case 'CCC': return 'P'
            case 'CCA': return 'P'
            case 'CCG': return 'P'
            case 'CAT': return 'H'
            case 'CAC': return 'H'
            case 'CAA': return 'Q'
            case 'CAG': return 'Q'
            case 'CGT': return 'R'
            case 'CGC': return 'R'
            case 'CGA': return 'R'
            case 'CGG': return 'R'
            case 'ATT': return 'I'
            case 'ATC': return 'I'
            case 'ATA': return 'I'
            case 'ATG': return 'M'
            case 'ACT': return 'T'
            case 'ACC': return 'T'
            case 'ACA': return 'T'
            case 'ACG': return 'T'
            case 'AAT': return 'N'
            case 'AAC': return 'N'
            case 'AAA': return 'K'
            case 'AAG': return 'K'
            case 'AGT': return 'S'
            case 'AGC': return 'S'
            case 'AGA': return 'R'
            case 'AGG': return 'R'
            case 'GTT': return 'V'
            case 'GTC': return 'V'
            case 'GTA': return 'V'
            case 'GTG': return 'V'
            case 'GCT': return 'A'
            case 'GCC': return 'A'
            case 'GCA': return 'A'
            case 'GCG': return 'A'
            case 'GAT': return 'D'
            case 'GAC': return 'D'
            case 'GAA': return 'E'
            case 'GAG': return 'E'
            case 'GGT': return 'G'
            case 'GGC': return 'G'
            case 'GGA': return 'G'
            case 'GGG': return 'G'
            default:
                if (codonUpper.contains("N") && codonUpper.length() == 3)
                    return "#"
                else
                    return '?'
        }
    }

    static String translate(String seq) {
        def aaSeq = ""
        def oof = seq.size() % 3
        if (oof > 0) {
            def mid = (int) (seq.size() / 2)
            seq = seq.substring(0, mid) + ("?" * (3 - oof)) + seq.substring(mid, seq.length())
        }

        def leftEnd = -1, rightEnd = -1
        for (int i = 0; i <= seq.size() - 3; i += 3) {
            def codon = seq.substring(i, i + 3)
            if (codon.contains("?")) {
                leftEnd = i
                break
            }
            aaSeq += codon2aa(codon)
        }

        if (oof == 0)
            return aaSeq

        def aaRight = ""
        for (int i = seq.size(); i >= 3; i -= 3) {
            def codon = seq.substring(i - 3, i)
            if (codon.contains("?")) {
                rightEnd = i
                break
            }
            aaRight += codon2aa(codon)
        }

        return aaSeq + seq.substring(leftEnd, rightEnd).toLowerCase() + aaRight.reverse()
    }


    @Override
    String toString() {
        [vSegment, dSegment, jSegment,
         cdr1start, cdr1end,
         cdr2start, cdr2end,
         cdr3start, cdr3end].join("\t")
    }
}
