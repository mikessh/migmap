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

package com.antigenomics.migmap.blast

import com.antigenomics.migmap.Util
import com.antigenomics.migmap.mapping.Cdr3Markup
import com.antigenomics.migmap.mapping.Truncations
import groovy.transform.CompileStatic

@CompileStatic
public class PSegmentSearcher {
    final static int P_SEGM_MIN_LEN = 3

    static PSegments search(Cdr3Markup cdr3Markup, Truncations truncations, String cdr3seq) {
        boolean dFound = cdr3Markup.dStart >= 0
        int vp3 = -1, dp5 = -1, dp3 = -1, jp5 = -1

        if (truncations.vDel == 0) {
            vp3 = searchRight(cdr3Markup.vEnd + 1, dFound ? cdr3Markup.dStart : cdr3Markup.jStart, 0, cdr3seq)
        }

        if (truncations.dDel5 == 0) {
            dp5 = searchLeft(cdr3Markup.dStart - 1, cdr3Markup.vEnd, cdr3Markup.dEnd, cdr3seq)
        }

        if (truncations.dDel3 == 0) {
            dp3 = searchRight(cdr3Markup.dEnd + 1, cdr3Markup.jStart, cdr3Markup.dStart, cdr3seq)
        }

        if (truncations.jDel == 0) {
            jp5 = searchLeft(cdr3Markup.jStart - 1, dFound ? cdr3Markup.vEnd : cdr3Markup.dEnd, cdr3seq.length(), cdr3seq)
        }

        new PSegments(
                vp3 - cdr3Markup.vEnd >= P_SEGM_MIN_LEN ? vp3 : -1,
                cdr3Markup.dStart - dp5 >= P_SEGM_MIN_LEN ? dp5 : -1,
                dp3 - cdr3Markup.dEnd >= P_SEGM_MIN_LEN ? dp3 : -1,
                cdr3Markup.jStart - jp5 >= P_SEGM_MIN_LEN ? jp5 : -1)
    }

    static int searchRight(int from, int to, int from2, String cdr3seq) {
        int lastMatch = -1, lastMm = -2
        for (int i = from; i < to; i++) {
            int j = 2 * from - i - 1
            if (j < from2)
                break
            if (cdr3seq.charAt(i) == Util.compl(cdr3seq.charAt(j))) {
                lastMatch = i
            } else {
                if (i - lastMm <= 2) { // no more than 2 consequent mismatches or 2 mismatches separated by 1 base
                    lastMatch = lastMm - 1
                    break
                }
                lastMm = i
            }
        }
        lastMatch
    }

    static int searchLeft(int from, int to, int from2, String cdr3seq) {
        int lastMatch = -1, lastMm = cdr3seq.length() + 1
        for (int i = from; i > to; i--) {
            int j = 2 * from - i + 1
            if (j >= from2)
                break
            if (cdr3seq.charAt(i) == Util.compl(cdr3seq.charAt(j))) {
                lastMatch = i
            } else {
                if (lastMm - i <= 2) { // no more than 2 consequent mismatches or 2 mismatches separated by 1 base
                    lastMatch = lastMm + 1
                    break
                }
                lastMm = i
            }
        }
        lastMatch
    }

}
