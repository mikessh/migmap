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
import groovy.transform.CompileStatic

@CompileStatic
public class PSegmentSearcher {
    final static int P_SEGM_MIN_LEN = 2

    static PSegments search(Cdr3Markup cdr3Markup, String cdr3seq) {
        int lastVMatch = -1, lastJMatch = -1, lastMm = -2

        for (int i = cdr3Markup.vEnd + 1;
             i < (cdr3Markup.dStart < 0 ? cdr3Markup.jStart : cdr3Markup.dStart); i++) {
            int j = 2 * cdr3Markup.vEnd - i + 1
            if (j < 0)
                break
            if (cdr3seq.charAt(i) == Util.compl(cdr3seq.charAt(j))) {
                lastVMatch = i
            } else {
                if (i - lastMm <= 2) { // no more than 2 consequent mismatches or 2 mismatches separated by 1 base
                    lastVMatch = lastMm - 1
                    break
                }
                lastMm = i
            }
        }

        lastMm = cdr3seq.length() + 1

        for (int i = cdr3Markup.jStart - 1;
             i > (cdr3Markup.dEnd < 0 ? cdr3Markup.dEnd : cdr3Markup.vEnd);
             i--) {
            int j = 2 * cdr3Markup.jStart - i - 1
            if (j >= cdr3seq.length())
                break
            if (cdr3seq.charAt(i) == Util.compl(cdr3seq.charAt(j))) {
                //println(cdr3seq.charAt(i).toString() + " " + Util.compl(cdr3seq.charAt(j)).toString() + "  ")
                lastJMatch = i
            } else {
                //println(cdr3seq.charAt(i).toString() + " " + Util.compl(cdr3seq.charAt(j)).toString() + " *")
                if (lastMm - i <= 2) { // no more than 2 consequent mismatches or 2 mismatches separated by 1 base
                    lastJMatch = lastMm + 1
                    break
                }
                lastMm = i
            }
        }

        new PSegments(lastVMatch - cdr3Markup.vEnd >= P_SEGM_MIN_LEN ? lastVMatch : -1,
                cdr3Markup.jStart - lastJMatch >= P_SEGM_MIN_LEN ? lastJMatch : -1)
    }
}
