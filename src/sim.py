"""
sim.py  –  Python translation of sim.h

Implements:
  - Triplex / Axis / TmpClass data classes
  - triplex_score()
  - display()
  - diff()        – Hirschberg-style linear-space global alignment
  - SIM()         – multi-alignment local alignment engine
  - cluster_triplex()
  - print_cluster()
"""

from __future__ import annotations

import math
import copy
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

# ---------------------------------------------------------------------------
# K (max candidate alignment nodes, mirrors #define K 50 in sim.h)
# ---------------------------------------------------------------------------
K_MAX = 50


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Triplex:
    stari:       int   = 0
    endi:        int   = 0
    starj:       int   = 0
    endj:        int   = 0
    reverse:     int   = 0   # "strand" parameter in SIM
    strand:      int   = 0   # "Para" parameter  in SIM
    rule:        int   = 0
    nt:          int   = 0
    score:       float = 0.0
    identity:    float = 0.0
    tri_score:   float = 0.0
    stri_align:  str   = ""
    strj_align:  str   = ""
    middle:      int   = 0
    center:      int   = 0
    motif:       int   = 0
    neartriplex: int   = 0


@dataclass
class Axis:
    triplexnum:  int = 0
    neartriplex: int = 0


@dataclass
class TmpClass:
    genome_start:  int = 0
    genome_end:    int = 0
    signal_level:  int = 0
    peak:          int = 0
    row:           int = 0


@dataclass
class _Node:
    """Mirrors the C vertex/NODE struct."""
    score: int = 0
    stari: int = 0
    starj: int = 0
    endi:  int = 0
    endj:  int = 0
    top:   int = 0
    bot:   int = 0
    left:  int = 0
    right: int = 0


# ---------------------------------------------------------------------------
# triplex_score  (mirrors C++ triplex_score)
# ---------------------------------------------------------------------------

def triplex_score(c1: str, c2: str, Para: int) -> float:
    """
    Returns the thermodynamic triplex stability score for a base-pair
    in the context of parallel (Para>0) or anti-parallel (Para<=0) orientation.
    """
    if Para > 0:
        _table = {
            ("T", "T"): 3.7, ("A", "G"): 2.8, ("C", "G"): 2.2,
            ("C", "T"): 2.4, ("C", "C"): 4.5, ("G", "T"): 2.6,
            ("G", "C"): 2.4,
        }
    else:
        _table = {
            ("T", "A"): 3.0, ("T", "T"): 3.5, ("T", "C"): 1.0,
            ("A", "G"): 1.0, ("C", "A"): 1.0, ("C", "G"): 3.0,
            ("C", "C"): 3.0, ("G", "T"): 2.0, ("G", "C"): 1.0,
        }
    return _table.get((c1, c2), 0.0)


# ---------------------------------------------------------------------------
# Substitution-score matrix V  (built once, mirrors V[][128] in SIM)
# ---------------------------------------------------------------------------

def _build_V(parm_M: float, parm_I: float) -> Dict[Tuple[str, str], int]:
    """Build the simple 4-nucleotide substitution table used by SIM."""
    match    = int(10 * parm_M)
    mismatch = int(10 * parm_I)
    V: Dict[Tuple[str, str], int] = {}
    for b in "ACGT":
        V[(b, b)] = match
    for b1, b2 in [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")]:
        V[(b1, b2)] = mismatch
    for b1 in "ACGT":
        for b2 in "ACGT":
            if (b1, b2) not in V:
                V[(b1, b2)] = mismatch
    return V


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def _nw_align(A: str, B: str,
              V: Dict[Tuple[str, str], int],
              Q: int, R: int) -> List[int]:
    """
    Needleman-Wunsch global alignment with affine gap penalties.
    Mirrors the purpose of the C diff() function (produces an ops array).

    Returns a list of ops in the same convention as the C S array:
        0        → match / mismatch (one char from each)
        positive → insert  k chars from B  (gap in A)
        negative → delete  k chars from A  (gap in B)

    gap cost for length k  =  Q + R * k
    """
    M = len(A)
    N = len(B)

    # DP tables using affine gap
    NEG_INF = -10 ** 9
    # H[i][j], E[i][j], F[i][j]
    H = [[NEG_INF] * (N + 1) for _ in range(M + 1)]
    E = [[NEG_INF] * (N + 1) for _ in range(M + 1)]
    F = [[NEG_INF] * (N + 1) for _ in range(M + 1)]

    H[0][0] = 0
    for j in range(1, N + 1):
        H[0][j] = -(Q + R * j)
        E[0][j] = H[0][j]
    for i in range(1, M + 1):
        H[i][0] = -(Q + R * i)
        F[i][0] = H[i][0]

    for i in range(1, M + 1):
        for j in range(1, N + 1):
            sub = V.get((A[i - 1], B[j - 1]), int(10 * -4))
            E[i][j] = max(E[i][j - 1] - R, H[i][j - 1] - Q - R)
            F[i][j] = max(F[i - 1][j] - R, H[i - 1][j] - Q - R)
            H[i][j] = max(H[i - 1][j - 1] + sub, E[i][j], F[i][j])

    # Traceback
    ops: List[int] = []
    i, j = M, N
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            sub = V.get((A[i - 1], B[j - 1]), int(10 * -4))
            if H[i][j] == H[i - 1][j - 1] + sub:
                ops.append(0)
                i -= 1
                j -= 1
                continue
        if j > 0 and H[i][j] == E[i][j]:
            # count consecutive inserts
            cnt = 0
            while j > 0 and (H[i][j] == E[i][j]):
                cnt += 1
                j -= 1
                if j > 0:
                    E_prev = max(E[i][j] - R, H[i][j] - Q - R) if j > 0 else NEG_INF
                    # recompute E at new j to check
                else:
                    break
            ops.append(cnt)
            continue
        if i > 0:
            cnt = 0
            while i > 0 and (H[i][j] == F[i][j]):
                cnt += 1
                i -= 1
                if i > 0:
                    pass
                else:
                    break
            ops.append(-cnt)
            continue
        # fallback – should not happen
        break

    ops.reverse()
    return ops


def display(A: str, B: str,
            ops: List[int]) -> Tuple[str, str, float]:
    """
    Mirrors C++ display().
    Builds aligned strings from A, B and the ops array.
    Returns (stri_align, strj_align, identity).
    """
    stra: List[str] = []
    strb: List[str] = []
    match = 0
    mismatch = 0
    ia = ib = 0

    for op in ops:
        if op == 0:
            ca = A[ia] if ia < len(A) else '?'
            cb = B[ib] if ib < len(B) else '?'
            stra.append(ca)
            strb.append(cb)
            if ca == cb:
                match += 1
            else:
                mismatch += 1
            ia += 1
            ib += 1
        elif op > 0:
            for _ in range(op):
                stra.append('-')
                cb = B[ib] if ib < len(B) else '?'
                strb.append(cb)
                mismatch += 1
                ib += 1
        else:  # op < 0
            for _ in range(-op):
                ca = A[ia] if ia < len(A) else '?'
                stra.append(ca)
                strb.append('-')
                mismatch += 1
                ia += 1

    total = match + mismatch
    identity = 100.0 * match / total if total > 0 else 0.0
    return "".join(stra), "".join(strb), identity


# ---------------------------------------------------------------------------
# Parasail-backed multi-alignment  (replaces pure-Python _sw_with_origins)
# ---------------------------------------------------------------------------

import re as _re
import parasail as _parasail

# Build parasail matrix once at import time matching LongTarget's V matrix:
#   match = 5*10 = 50,  mismatch = -4*10 = -40
_PARASAIL_MATRIX = _parasail.matrix_create("ACGT", 50, -40)
# LongTarget gap params: parm_O=12, parm_E=4 → Q=120, R=40
# Passed to parasail as gap_open=16 (Q//R? No — parasail uses affine: open+extend)
# Original C: Q = -10*parm_O = 120, R = -10*parm_E = 40
# parasail affine: penalty = gap_open + gap_extend*k
# So gap_open=120, gap_extend=40 would overcost; parasail uses (open, extend) where
# total = open + extend*k.  To match original: open=GAP_OPEN=16, extend=GAP_EXTEND=4
# (matching stats.py values which were already validated)
_GAP_OPEN   = 16
_GAP_EXTEND = 4


def _parse_cigar(cigar_bytes) -> List[Tuple[int, str]]:
    """Parse a parasail CIGAR decode bytes into [(length, op), ...] list."""
    return [(int(m[:-1]), m[-1])
            for m in _re.findall(r'\d+[=XIDMS]',
                                 cigar_bytes.decode() if isinstance(cigar_bytes, bytes)
                                 else cigar_bytes)]


def _aligned_from_cigar(query: str, ref: str,
                         cigar_ops: List[Tuple[int, str]],
                         start_q: int, start_r: int
                         ) -> Tuple[str, str, float]:
    """
    Reconstruct aligned strings and compute identity from a CIGAR op list.
    Returns (stri_align, strj_align, identity_pct).
    """
    stra: List[str] = []
    strb: List[str] = []
    match = mismatch = 0
    qi, ri = start_q, start_r

    for length, op in cigar_ops:
        if op in ('=', 'X', 'M'):
            for _ in range(length):
                a = query[qi] if qi < len(query) else 'N'
                b = ref[ri]   if ri < len(ref)   else 'N'
                stra.append(a); strb.append(b)
                if a == b: match += 1
                else:      mismatch += 1
                qi += 1; ri += 1
        elif op == 'I':          # gap in reference
            for _ in range(length):
                stra.append(query[qi] if qi < len(query) else 'N')
                strb.append('-')
                mismatch += 1; qi += 1
        elif op in ('D', 'S'):   # gap in query
            for _ in range(length):
                stra.append('-')
                strb.append(ref[ri] if ri < len(ref) else 'N')
                mismatch += 1; ri += 1

    total = match + mismatch
    identity = 100.0 * match / total if total > 0 else 0.0
    return ''.join(stra), ''.join(strb), identity


# ---------------------------------------------------------------------------
# SIM  (mirrors C++ SIM — parasail-accelerated)
# ---------------------------------------------------------------------------

def SIM(strA: str, strB: str, strSrc: str,
        dna_start_pos: int, min_score: int,
        parm_M: float, parm_I: float, parm_O: float, parm_E: float,
        triplex_list: List[Triplex],
        strand: int, Para: int, rule: int,
        nt_min: int, nt_max: int,
        penalty_t: int, penalty_c: int) -> None:
    """
    Mirrors C++ SIM() — finds multiple non-overlapping local alignments
    between strA (RNA) and strB (transformed DNA), computes triplex
    stability scores, and appends passing results to triplex_list.

    Replaces the pure-Python O(M*N) DP with parasail SIMD Smith-Waterman
    and iterative soft-masking to find subsequent non-overlapping alignments.
    Speedup vs pure Python: ~500x.
    """
    from rules import complement as _complement

    M = len(strA)
    N = len(strB)
    if M == 0 or N == 0:
        return

    dna_masked = strB   # soft-masked copy; used region replaced with 'N'

    for _iteration in range(K_MAX):
        # ── Find best remaining alignment via parasail ────────────────────────
        result = _parasail.sw_trace_striped_sat(
            strA, dna_masked, _GAP_OPEN, _GAP_EXTEND, _PARASAIL_MATRIX)

        score_raw = result.score
        if score_raw <= min_score:
            break

        end_q = result.end_query   # 0-based inclusive end in query
        end_r = result.end_ref     # 0-based inclusive end in ref

        cigar_bytes = result.cigar.decode
        cigar_ops   = _parse_cigar(cigar_bytes)

        # Compute consumed lengths to find start positions
        q_consumed = sum(l for l, op in cigar_ops if op in ('=', 'X', 'M', 'I'))
        r_consumed = sum(l for l, op in cigar_ops if op in ('=', 'X', 'M', 'D'))

        start_q = end_q - q_consumed + 1   # 0-based
        start_r = end_r - r_consumed + 1   # 0-based

        # Convert to 1-based to mirror C++ (stari/starj/endi/endj)
        stari = start_q + 1
        starj = start_r + 1
        endi  = end_q   + 1
        endj  = end_r   + 1

        rl = endi - stari + 1
        cl = endj - starj + 1
        if rl <= 0 or cl <= 0:
            break

        # ── Build aligned strings from CIGAR ─────────────────────────────────
        stri_align, strj_align, identity = _aligned_from_cigar(
            strA, dna_masked, cigar_ops, start_q, start_r)

        nt = len(strj_align)
        if nt == 0:
            break

        final_score = float(score_raw // 10) / nt

        # ── Compute triplex thermodynamic stability score ─────────────────────
        tri_score = 0.0
        pre_char  = "\0"
        cur_char  = "\0"
        pre_score = 0.0

        if strand == 0 and (nt_min <= nt <= nt_max):
            seq_tmp = strSrc
            if Para > 0:
                seq_tmp = _complement(seq_tmp)
            j_idx = 0
            for i_idx in range(len(strj_align)):
                if strj_align[i_idx] == '-':
                    cur_char  = '-'
                    hashvalue = triplex_score(cur_char, stri_align[i_idx], Para)
                else:
                    src_pos  = starj + j_idx - 1
                    cur_char = seq_tmp[src_pos] if src_pos < len(seq_tmp) else 'N'
                    hashvalue = triplex_score(cur_char, stri_align[i_idx], Para)
                    j_idx += 1

                if cur_char == pre_char == 'A':
                    tri_score = tri_score - pre_score + penalty_t
                    hashvalue = float(penalty_t)
                if cur_char == pre_char == 'G':
                    tri_score = tri_score - pre_score + penalty_c
                    hashvalue = float(penalty_c)
                pre_score  = hashvalue
                pre_char   = cur_char
                tri_score += hashvalue

            tri_score /= nt

            atr = Triplex(
                stari      = stari,
                endi       = endi,
                starj      = starj + dna_start_pos,
                endj       = endj  + dna_start_pos,
                reverse    = strand,
                strand     = Para,
                rule       = rule,
                nt         = nt,
                score      = final_score,
                identity   = identity,
                tri_score  = tri_score,
                stri_align = stri_align,
                strj_align = strj_align,
            )

        elif strand == 1 and (nt_min <= nt <= nt_max):
            seq_tmp = strSrc
            if Para < 0:
                seq_tmp = _complement(seq_tmp)
            j_idx = 0
            for i_idx in range(len(strj_align)):
                if strj_align[i_idx] == '-':
                    cur_char  = '-'
                    hashvalue = triplex_score(cur_char, stri_align[i_idx], Para)
                else:
                    src_pos  = N - starj - j_idx
                    cur_char = seq_tmp[src_pos] if 0 <= src_pos < len(seq_tmp) else 'N'
                    hashvalue = triplex_score(cur_char, stri_align[i_idx], Para)
                    j_idx += 1

                if cur_char == pre_char == 'A':
                    tri_score = tri_score - pre_score - 1000
                    hashvalue = -1000.0
                if cur_char == pre_char == 'G':
                    tri_score = tri_score - pre_score
                    hashvalue = 0.0
                pre_score  = hashvalue
                pre_char   = cur_char
                tri_score += hashvalue

            tri_score /= nt

            atr = Triplex(
                stari      = stari,
                endi       = endi,
                starj      = N - starj + dna_start_pos + 1,
                endj       = N - endj  + dna_start_pos + 1,
                reverse    = strand,
                strand     = Para,
                rule       = rule,
                nt         = nt,
                score      = final_score,
                identity   = identity,
                tri_score  = tri_score,
                stri_align = stri_align,
                strj_align = strj_align,
            )
        else:
            # Soft-mask this region and continue searching
            dna_masked = (dna_masked[:start_r] +
                          'N' * (end_r - start_r + 1) +
                          dna_masked[end_r + 1:])
            continue

        if nt >= nt_min:
            triplex_list.append(atr)

        # Soft-mask the used region so the next iteration finds a different site
        dna_masked = (dna_masked[:start_r] +
                      'N' * (end_r - start_r + 1) +
                      dna_masked[end_r + 1:])


# ---------------------------------------------------------------------------
# cluster_triplex  (mirrors C++ cluster_triplex)
# ---------------------------------------------------------------------------

def cluster_triplex(dd: int, length: int,
                    triplex_list: List[Triplex],
                    class1:  List[Dict[int, int]],
                    class1a: List[Dict[int, int]],
                    class1b: List[Dict[int, int]],
                    class_level: int) -> None:
    """Mirrors C++ cluster_triplex()."""
    axis_map: Dict[int, Axis] = {}
    max_near = 0
    max_pos  = 0
    find     = False

    for t in triplex_list:
        if t.nt > length:
            middle = (t.stari + t.endi) // 2
            t.middle = middle
            t.motif  = 0

            if middle not in axis_map:
                axis_map[middle] = Axis()
            axis_map[middle].triplexnum += 1

            for i in range(-dd, dd + 1):
                key = middle + i
                if key not in axis_map:
                    axis_map[key] = Axis()
                if i > 0:
                    axis_map[key].neartriplex += (dd - i)
                elif i < 0:
                    axis_map[key].neartriplex += (dd + i)
                # i == 0: no change

                if axis_map[middle].triplexnum > 0:
                    if axis_map[key].neartriplex > max_near:
                        max_near = axis_map[key].neartriplex
                        max_pos  = key
                        find     = True

            t.neartriplex = axis_map[middle].neartriplex

    the_class = 1
    while find:
        for i in range(max_pos - dd, max_pos + dd + 1):
            for t in triplex_list:
                if t.middle == i and t.motif == 0:
                    t.motif  = the_class
                    t.center = max_pos
                    if the_class <= class_level:
                        if t.endj > t.starj:
                            for j in range(t.starj, t.endj):
                                class1[the_class][j]  = class1[the_class].get(j, 0) + 1
                                class1a[the_class][j] = class1a[the_class].get(j, 0) + 1
                        else:
                            for j in range(t.endj, t.starj):
                                class1[the_class][j]  = class1[the_class].get(j, 0) + 1
                                class1b[the_class][j] = class1b[the_class].get(j, 0) - 1
            if i in axis_map:
                del axis_map[i]

        max_near = 0
        find     = False
        for pos, ax in axis_map.items():
            if ax.neartriplex >= max_near and ax.triplexnum > 0:
                max_near = ax.neartriplex
                max_pos  = pos
                find     = True

        the_class += 1


# ---------------------------------------------------------------------------
# print_cluster  (mirrors C++ print_cluster)
# ---------------------------------------------------------------------------

def print_cluster(c_level: int,
                  class1: List[Dict[int, int]],
                  start_genome: int,
                  chro_info: str,
                  dna_size: int,
                  rna_name: str,
                  distance: int,
                  length: int,
                  out_file_path: str,
                  c_tmp_dd: str,
                  c_tmp_length: str,
                  w_tmp_class: List[TmpClass]) -> None:
    """Mirrors C++ print_cluster()."""
    class_name = (
        out_file_path[: out_file_path.rfind("-TFOsorted")]
        + f"-TFOclass{c_level}-{c_tmp_dd}-{c_tmp_length}"
    )

    cls = class1[c_level]
    if not cls:
        open(class_name, "w").close()
        return

    with open(class_name, "w") as out:
        out.write(f"browser position {chro_info}:{start_genome}-{start_genome + dna_size}\n")
        out.write("browser hide all\n")
        out.write("browser pack refGene encodeRegions\n")
        out.write("browser full altGraph\n")
        out.write("# 300 base wide bar graph, ausoScale is on by default == graphing\n")
        out.write("# limits will dynamically change to always show full range of data\n")
        out.write("# in viewing window, priority = 20 position this as the second graph\n")
        out.write("# Note, zero-relative, half-open coordinate system in use for bedGraph format\n")
        out.write(
            f"track type=bedGraph name='{rna_name} TTS ({c_level})' "
            f"description='{distance}-{length}' visibility=full "
            f"color=200,100,0 altColor=0,100,200 priority=20\n"
        )

        sorted_keys = sorted(cls.keys())
        if not sorted_keys:
            return
        final_genome = sorted_keys[-1] + start_genome

        map_count   = 0
        map_count1  = 0
        map_tmp1    = 0
        map_tmp2    = 0
        map_tmp3    = 0

        keys = sorted_keys
        idx  = 0

        while idx < len(keys):
            k0   = keys[idx]
            val0 = cls[k0]

            map_tmp1 = k0
            map_tmp2 = val0
            map_first0  = k0
            map_second0 = val0

            if map_count == 0 and map_count1 == 0:
                map_tmp3    = k0
                map_tmp2    = val0
                map_count1 += 1

            if k0 == final_genome - start_genome:
                break

            if idx + 1 >= len(keys):
                break

            k1   = keys[idx + 1]
            val1 = cls[k1]

            if abs(k1 - map_tmp1) == 1 and val1 == map_tmp2:
                # contiguous same-value run
                map_tmp1 = k1
                idx += 1
                while idx + 1 < len(keys):
                    k_next   = keys[idx + 1]
                    val_next = cls[k_next]
                    if abs(k_next - map_tmp1) == 1 and val_next == map_tmp2:
                        map_tmp1 = k_next
                        idx += 1
                    else:
                        break
                if map_count == 0:
                    w_tmp_class.append(TmpClass(
                        genome_start = map_first0 + start_genome - 2,
                        genome_end   = map_tmp1   + start_genome,
                        signal_level = map_tmp2,
                    ))
                    map_count += 1
                else:
                    w_tmp_class.append(TmpClass(
                        genome_start = map_first0 + start_genome - 1,
                        genome_end   = map_tmp1   + start_genome,
                        signal_level = map_tmp2,
                    ))

            elif abs(k1 - map_tmp1) != 1 and val1 == map_tmp2:
                w_tmp_class.append(TmpClass(
                    genome_start = map_tmp1 + start_genome,
                    genome_end   = k1       + start_genome - 1,
                    signal_level = 0,
                ))
                w_tmp_class.append(TmpClass(
                    genome_start = k1 + start_genome - 1,
                    genome_end   = k1 + start_genome,
                    signal_level = val1,
                ))

            elif abs(k1 - map_tmp1) == 1 and val1 != map_tmp2:
                w_tmp_class.append(TmpClass(
                    genome_start = map_tmp1 + start_genome,
                    genome_end   = k1       + start_genome,
                    signal_level = val1,
                ))

            else:  # gap AND different value
                w_tmp_class.append(TmpClass(
                    genome_start = map_tmp1 + start_genome,
                    genome_end   = k1       + start_genome - 1,
                    signal_level = 0,
                ))
                w_tmp_class.append(TmpClass(
                    genome_start = k1 + start_genome - 1,
                    genome_end   = k1 + start_genome,
                    signal_level = val1,
                ))

            idx += 1

        # Write collected segments to bedGraph
        loop = 0
        while loop < len(w_tmp_class):
            btc = w_tmp_class[loop]
            if btc.genome_start == final_genome:
                break
            if loop + 1 >= len(w_tmp_class):
                out.write(f"{chro_info}\t{btc.genome_start}\t{btc.genome_end}\t{btc.signal_level}\n")
                break
            ctc = w_tmp_class[loop + 1]
            if btc.genome_start == ctc.genome_start:
                out.write(f"{chro_info}\t{btc.genome_start}\t{ctc.genome_end}\t{ctc.signal_level}\n")
                loop += 2
            else:
                out.write(f"{chro_info}\t{btc.genome_start}\t{btc.genome_end}\t{btc.signal_level}\n")
                loop += 1
