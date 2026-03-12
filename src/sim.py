"""
sim.py – Python translation of sim.h
Triplex DNA/RNA local alignment (SIM algorithm)
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Dict, Set, Optional
import copy
from rules import complement

# ─────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────

K = 50  # Maximum number of alignment nodes kept in LIST


# ─────────────────────────────────────────────────────────────
# Data-structures  (C++ structs → Python dataclasses)
# ─────────────────────────────────────────────────────────────

@dataclass
class Triplex:
    stari:      int   = 0
    endi:       int   = 0
    starj:      int   = 0
    endj:       int   = 0
    reverse:    int   = 0
    strand:     int   = 0
    rule:       int   = 0
    nt:         int   = 0
    score:      float = 0.0
    identity:   float = 0.0
    tri_score:  float = 0.0
    stri_align: str   = ""
    strj_align: str   = ""
    middle:     int   = 0
    center:     int   = 0
    motif:      int   = 0
    neartriplex:int   = 0


@dataclass
class Node:
    """Corresponds to the NODE / vertex struct."""
    SCORE: int = 0
    STARI: int = 0
    STARJ: int = 0
    ENDI:  int = 0
    ENDJ:  int = 0
    TOP:   int = 0
    BOT:   int = 0
    LEFT:  int = 0
    RIGHT: int = 0


@dataclass
class Axis:
    triplexnum:  int = 0
    neartriplex: int = 0


@dataclass
class TmpClass:
    genome_start: int = 0
    genome_end:   int = 0
    signal_level: int = 0
    peak:         int = 0
    row:          int = 0


# ─────────────────────────────────────────────────────────────
# Utility
# ─────────────────────────────────────────────────────────────

def gap_cost(k: int, Q: int, R: int) -> int:
    """Affine gap cost: g(k) = Q + R*k  (0 for k<=0)."""
    return 0 if k <= 0 else Q + R * k


def order(ss1: int, xx1: int, yy1: int,
          ss2: int, xx2: int, yy2: int):
    """
    C++ ORDER macro:  if ss1 < ss2 replace (ss1,xx1,yy1) with (ss2,xx2,yy2);
    break ties by maximising xx then yy.
    Returns (ss1, xx1, yy1) after possible update.
    """
    if ss1 < ss2:
        return ss2, xx2, yy2
    if ss1 == ss2:
        if xx1 < xx2:
            return ss1, xx2, yy2
        if xx1 == xx2 and yy1 < yy2:
            return ss1, xx1, yy2
    return ss1, xx1, yy1


# ─────────────────────────────────────────────────────────────
# triplex_score
# ─────────────────────────────────────────────────────────────

def triplex_score(c1: str, c2: str, Para: int) -> float:
    """Score a base-pair in a triplex according to the orientation (Para)."""
    if Para > 0:
        table = {
            ('T', 'T'): 3.7, ('A', 'G'): 2.8, ('C', 'G'): 2.2,
            ('C', 'T'): 2.4, ('C', 'C'): 4.5, ('G', 'T'): 2.6,
            ('G', 'C'): 2.4,
        }
    else:
        table = {
            ('T', 'A'): 3.0, ('T', 'T'): 3.5, ('T', 'C'): 1.0,
            ('A', 'G'): 1.0, ('C', 'A'): 1.0, ('C', 'G'): 3.0,
            ('C', 'C'): 3.0, ('G', 'T'): 2.0, ('G', 'C'): 1.0,
        }
    return table.get((c1, c2), 0.0)


# ─────────────────────────────────────────────────────────────
# addnode
# ─────────────────────────────────────────────────────────────

def addnode(c: int, ci: int, cj: int, i: int, j: int,
            LIST: List[Node], numnode_ref: List[int]) -> int:
    """
    Insert or update a node in LIST.
    numnode_ref is a single-element list so it can be mutated (simulates C pointer).
    Returns 1.
    """
    numnode = numnode_ref[0]
    found = False
    most = 0

    for d in range(numnode):
        most = d
        if LIST[most].STARI == ci and LIST[most].STARJ == cj:
            found = True
            break

    if found:
        if LIST[most].SCORE < c:
            LIST[most].SCORE = c
            LIST[most].ENDI  = i
            LIST[most].ENDJ  = j
        if LIST[most].TOP   > i: LIST[most].TOP   = i
        if LIST[most].BOT   < i: LIST[most].BOT   = i
        if LIST[most].LEFT  > j: LIST[most].LEFT  = j
        if LIST[most].RIGHT < j: LIST[most].RIGHT = j
    else:
        low = 0
        if numnode == K:
            for d in range(1, numnode):
                if LIST[d].SCORE < LIST[low].SCORE:
                    low = d
            most = low
        else:
            most = numnode
            numnode_ref[0] += 1

        LIST[most].SCORE = c
        LIST[most].STARI = ci
        LIST[most].STARJ = cj
        LIST[most].ENDI  = i
        LIST[most].ENDJ  = j
        LIST[most].TOP   = i
        LIST[most].BOT   = i
        LIST[most].LEFT  = j
        LIST[most].RIGHT = j

    return 1


# ─────────────────────────────────────────────────────────────
# no_cross
# ─────────────────────────────────────────────────────────────

def no_cross(LIST: List[Node], numnode: int,
             m1: int, mm: int, n1: int, nn: int,
             prl: List[int], pcl: List[int]) -> int:
    """
    Check whether the current alignment region overlaps any node that would
    cause a crossing.  Returns 1 (no crossing) or 0 (crossing found).
    prl / pcl are single-element lists (simulate C long* pointers).
    """
    for i in range(numnode):
        cur = LIST[i]
        if (cur.STARI <= mm and cur.STARJ <= nn and
                cur.BOT  >= m1 - 1 and cur.RIGHT >= n1 - 1 and
                (cur.STARI < prl[0] or cur.STARJ < pcl[0])):
            if cur.STARI < prl[0]: prl[0] = cur.STARI
            if cur.STARJ < pcl[0]: pcl[0] = cur.STARJ
            break
    else:
        return 1   # loop completed without break → no crossing
    return 0


# ─────────────────────────────────────────────────────────────
# diff  –  divide-and-conquer alignment (recursive)
# ─────────────────────────────────────────────────────────────

class _DiffState:
    """
    Bundles the mutable state that travels through every recursive diff call:
      S      – the edit-script array
      sapp   – write cursor into S  (replaces the C 'sapp' pointer)
      plast  – last value written   (replaces C '*plast')
    """
    __slots__ = ("S", "sapp", "plast")

    def __init__(self, size: int) -> None:
        self.S:     List[int] = [0] * size
        self.sapp:  int       = 0
        self.plast: int       = 0


def _del(k: int, pI: List[int], st: _DiffState) -> None:
    """C macro DEL(k) – record k deletions in the edit script."""
    pI[0] += k
    if st.plast < 0:
        st.S[st.sapp - 1] -= k
        st.plast = st.S[st.sapp - 1]
    else:
        st.S[st.sapp] = -k
        st.plast = -k
        st.sapp += 1


def _ins(k: int, pJ: List[int], st: _DiffState) -> None:
    """C macro INS(k) – record k insertions in the edit script."""
    pJ[0] += k
    if st.plast < 0:
        old = st.plast
        st.S[st.sapp - 1] = k
        st.S[st.sapp]     = old
        st.sapp += 1
    else:
        st.S[st.sapp] = k
        st.plast = k
        st.sapp += 1


def _rep(st: _DiffState) -> None:
    """C macro REP – record a replacement in the edit script."""
    st.S[st.sapp] = 0
    st.plast = 0
    st.sapp += 1


def diff(A: str, A_off: int,
         B: str, B_off: int,
         M: int, N: int,
         pI: List[int], pJ: List[int],
         tb: int, te: int, Q: int, R: int,
         st: _DiffState,
         V: List[List[int]],
         row: List[Set[int]],
         CC: List[int], DD: List[int],
         RR: List[int], SS: List[int]) -> int:
    """
    Hirschberg-style divide-and-conquer alignment.

    A[A_off+1 .. A_off+M]  and  B[B_off+1 .. B_off+N]  are aligned
    (1-based convention carried over from C, where the strings were
    prepended with a space character).

    pI, pJ  – mutable [int] wrappers for the running row/column counters.
    st      – shared edit-script state.
    V       – 128×128 substitution score matrix.
    row     – list of sets: row[i] holds j-values already used on row i.
    CC/DD/RR/SS – work arrays (pre-allocated, length ≥ N+1).
    """

    # ---- helpers to fetch characters (1-based) -------------------------
    def Ac(i: int) -> int: return ord(A[A_off + i])
    def Bc(j: int) -> int: return ord(B[B_off + j])

    # ---- base cases ----------------------------------------------------
    if N <= 0:
        if M > 0:
            _del(M, pI, st)
        return -gap_cost(M, Q, R)

    if M <= 1:
        if M <= 0:
            _ins(N, pJ, st)
            return -gap_cost(N, Q, R)

        if tb > te:
            tb = te
        midc = -(tb + R + gap_cost(N, Q, R))
        midj = 0
        va = V[Ac(1)]
        for j in range(1, N + 1):
            ii = pI[0] + 1
            jj = j + pJ[0]
            if jj not in row[ii]:                       # DIAG check
                c = va[Bc(j)] - (gap_cost(j - 1, Q, R) + gap_cost(N - j, Q, R))
                if c > midc:
                    midc = c
                    midj = j

        if midj == 0:
            _ins(N, pJ, st)
            _del(1, pI, st)
        else:
            if midj > 1:
                _ins(midj - 1, pJ, st)
            _rep(st)
            pI[0] += 1
            pJ[0] += 1
            row[pI[0]].add(pJ[0])
            if midj < N:
                _ins(N - midj, pJ, st)
        return midc

    # ---- divide --------------------------------------------------------
    midi = M // 2

    # Forward half: rows 0..midi
    CC[0] = 0
    t = -Q
    for j in range(1, N + 1):
        t -= R
        CC[j] = t
        DD[j] = t - Q

    t = -tb
    for i in range(1, midi + 1):
        s      = CC[0]
        t     -= R
        c      = t
        CC[0]  = c
        e      = c - Q
        va     = V[Ac(i)]
        for j in range(1, N + 1):
            # gap extend horizontally
            c_h = c - Q - R
            e_h = e - R
            e = c_h if c_h > e_h else e_h
            # gap extend vertically
            c_v = CC[j] - Q - R
            d   = DD[j] - R
            if c_v > d: d = c_v
            c = c_v
            # diagonal (substitution)
            ii = i + pI[0];  jj = j + pJ[0]
            if jj not in row[ii]:
                c = s + va[Bc(j)]
            # best of three
            if c < d: c = d
            if c < e: c = e
            s     = CC[j]
            CC[j] = c
            DD[j] = d

    DD[0] = CC[0]

    # Backward half: rows midi..M-1
    RR[N] = 0
    t = -Q
    for j in range(N - 1, -1, -1):
        t -= R
        RR[j] = t
        SS[j] = t - Q

    t = -te
    for i in range(M - 1, midi - 1, -1):
        s      = RR[N]
        t     -= R
        c      = t
        RR[N]  = c
        e      = c - Q
        va     = V[Ac(i + 1)]
        for j in range(N - 1, -1, -1):
            c_h = c - Q - R
            e_h = e - R
            e = c_h if c_h > e_h else e_h
            c_v = RR[j] - Q - R
            d   = SS[j] - R
            if c_v > d: d = c_v
            c = c_v
            ii = i + 1 + pI[0];  jj = j + 1 + pJ[0]
            if jj not in row[ii]:
                c = s + va[Bc(j + 1)]
            if c < d: c = d
            if c < e: c = e
            s     = RR[j]
            RR[j] = c
            SS[j] = d

    SS[N] = RR[N]

    # Find best midpoint column
    midc = CC[0] + RR[0]
    midj = 0
    typ  = 1
    for j in range(N + 1):
        c = CC[j] + RR[j]
        if c >= midc:
            if c > midc or (CC[j] != DD[j] and RR[j] == SS[j]):
                midc = c
                midj = j
    for j in range(N, -1, -1):
        c = DD[j] + SS[j] + Q
        if c > midc:
            midc = c
            midj = j
            typ  = 2

    # Conquer
    if typ == 1:
        diff(A, A_off,        B, B_off,        midi,     midj,     pI, pJ, tb, Q,  Q, R, st, V, row, CC, DD, RR, SS)
        diff(A, A_off + midi, B, B_off + midj, M - midi, N - midj, pI, pJ, Q,  te, Q, R, st, V, row, CC, DD, RR, SS)
    else:
        diff(A, A_off,            B, B_off,        midi - 1,     midj, pI, pJ, tb, 0,  Q, R, st, V, row, CC, DD, RR, SS)
        _del(2, pI, st)
        diff(A, A_off + midi + 1, B, B_off + midj, M - midi - 1, N - midj, pI, pJ, 0, te, Q, R, st, V, row, CC, DD, RR, SS)

    return midc


# ─────────────────────────────────────────────────────────────
# display  –  build alignment strings and compute identity
# ─────────────────────────────────────────────────────────────

def display(A: str, A_off: int,
            B: str, B_off: int,
            M: int, N: int,
            S: List[int]) -> tuple[str, str, float]:
    """
    Walk the edit script S and build the two aligned strings.
    Returns (stri_align, strj_align, identity_percent).

    A[A_off+1..A_off+M], B[B_off+1..B_off+N]  (1-based, matching C convention).
    """
    stra:       List[str] = []
    strb:       List[str] = []
    match       = 0
    mis_match   = 0
    i = j = 0
    s_idx = 0          # read cursor into S

    while i < M or j < N:
        # consume zero-ops (matches / mismatches)
        while i < M and j < N and S[s_idx] == 0:
            i += 1;  j += 1
            a_ch = A[A_off + i]
            b_ch = B[B_off + j]
            if a_ch == b_ch:
                match += 1
            else:
                mis_match += 1
            stra.append(a_ch)
            strb.append(b_ch)
            s_idx += 1

        if i < M or j < N:
            op = S[s_idx];  s_idx += 1
            if op > 0:                  # insertion into A
                for _ in range(op):
                    stra.append('-')
                    j += 1
                    strb.append(B[B_off + j])
                    mis_match += 1
            else:                       # deletion from B
                for _ in range(-op):
                    strb.append('-')
                    i += 1
                    stra.append(A[A_off + i])
                    mis_match += 1

    total    = match + mis_match
    identity = 100.0 * match / total if total else 0.0
    return ''.join(stra), ''.join(strb), identity


# ─────────────────────────────────────────────────────────────
# SIM  –  main local-alignment driver
# ─────────────────────────────────────────────────────────────

def SIM(strA: str, strB: str, strSrc: str,
        dnaStartPos: int, min_score: int,
        parm_M: float, parm_I: float, parm_O: float, parm_E: float,
        triplex_list: List[Triplex],
        strand: int, Para: int, rule: int,
        ntMin: int, ntMax: int,
        penaltyT: int, penaltyC: int) -> None:
    """
    Local alignment of strA (RNA / TFO) against strB (DNA duplex region).
    Appends Triplex records to triplex_list.
    """

    # ---- scoring matrix V[128][128] ------------------------------------
    V: List[List[int]] = [[0] * 128 for _ in range(128)]

    def set_v(c1: str, c2: str, val: int) -> None:
        V[ord(c1)][ord(c2)] = val

    match_val = int(10 * parm_M)
    trans_val = int(10 * parm_I)   # transitions
    tvers_val = int(10 * parm_I)   # transversions (same weight here)

    set_v('A','A', match_val);  set_v('C','C', match_val)
    set_v('G','G', match_val);  set_v('T','T', match_val)

    for pair in [('A','G'),('G','A'),('C','T'),('T','C')]:
        set_v(pair[0], pair[1], trans_val)
    for pair in [('A','C'),('A','T'),('C','A'),('C','G'),
                 ('G','C'),('G','T'),('T','A'),('T','G')]:
        set_v(pair[0], pair[1], tvers_val)

    Q: int = int(-10 * parm_O)    # gap-open penalty (positive)
    R: int = int(-10 * parm_E)    # gap-extend penalty (positive)

    # ---- prepend a space so strings become 1-indexed -------------------
    A = ' ' + strA          # A[1..M]
    B = ' ' + strB          # B[1..N]
    M = len(strA)
    N = len(strB)

    # ---- work arrays ---------------------------------------------------
    CC = [0] * (N + 1);  DD = [0] * (N + 1)
    RR = [0] * (N + 1);  SS = [0] * (N + 1)
    EE = [0] * (N + 1);  FF = [0] * (N + 1)
    HH = [0] * (M + 1);  WW = [0] * (M + 1)
    II = [0] * (M + 1);  JJ = [0] * (M + 1)
    XX = [0] * (M + 1);  YY = [0] * (M + 1)
    S_buf = [0] * (N + M + 2)     # edit-script buffer

    # row[i] holds j-values already used on row i (to avoid repeats)
    row: List[Set[int]] = [set() for _ in range(M + 2)]

    LIST:    List[Node] = [Node() for _ in range(K)]
    numnode: List[int]  = [0]

    # ─────────────── first pass: fill the DP table ──────────────────────
    for j in range(1, N + 1):
        CC[j] = 0
        RR[j] = 0
        EE[j] = j
        DD[j] = -(Q)
        SS[j] = 0
        FF[j] = j

    for i in range(1, M + 1):
        c = 0
        f = -(Q)
        ci = fi = i
        va = V[ord(A[i])]
        p  = 0
        pi = i - 1
        cj = fj = pj = 0

        for j in range(1, N + 1):
            f -= R
            c2 = c - Q - R
            f, fi, fj = order(f, fi, fj, c2, ci, cj)

            c  = CC[j] - Q - R
            ci = RR[j];  cj = EE[j]
            d  = DD[j] - R
            di = SS[j];  dj = FF[j]
            d, di, dj = order(d, di, dj, c, ci, cj)

            # DIAG
            c = 0
            if j not in row[i]:
                c = p + va[ord(B[j])]

            if c <= 0:
                c = 0;  ci = i;  cj = j
            else:
                ci = pi;  cj = pj

            c, ci, cj = order(c, ci, cj, d, di, dj)
            c, ci, cj = order(c, ci, cj, f, fi, fj)

            p  = CC[j];  CC[j] = c
            pi = RR[j];  pj = EE[j]
            RR[j] = ci;  EE[j] = cj
            DD[j] = d
            SS[j] = di;  FF[j] = dj

            if c > min_score:
                addnode(c, ci, cj, i, j, LIST, numnode)

    # ─────────────── process nodes (highest score first) ─────────────────
    for count in range(numnode[0] - 1, -1, -1):

        # find highest-score node
        best = 0
        for idx in range(1, numnode[0]):
            if LIST[idx].SCORE > LIST[best].SCORE:
                best = idx

        cur = copy.copy(LIST[best])

        # swap best to end and shrink list
        numnode[0] -= 1
        if best != numnode[0]:
            LIST[best] = copy.copy(LIST[numnode[0]])
            LIST[numnode[0]] = cur

        score = cur.SCORE
        stari = cur.STARI + 1
        starj = cur.STARJ + 1
        endi  = cur.ENDI
        endj  = cur.ENDJ
        m1 = cur.TOP;  mm = cur.BOT
        n1 = cur.LEFT; nn = cur.RIGHT

        rl = endi - stari + 1
        cl = endj - starj + 1
        pI = [stari - 1]
        pJ = [starj - 1]

        st = _DiffState(N + M + 2)
        diff(A, stari - 1, B, starj - 1, rl, cl,
             pI, pJ, Q, Q, Q, R,
             st, V, row, CC, DD, RR, SS)

        if score / 10.0 <= min_score:
            break

        nt = endi - stari + 1
        stri_align, strj_align, identity = display(
            A, stari - 1, B, starj - 1, rl, cl, st.S)
        nt = len(strj_align)

        atriplex  = Triplex()
        tri_score = 0.0
        final_score = 0.0

        if strand == 0 and ntMin <= nt <= ntMax:
            seqtmp = strSrc
            if Para > 0:
                seqtmp = complement(seqtmp)
            j_idx = 0
            hashvalue = prescore = 0.0
            prechar = curchar = ''
            for idx_i in range(len(strj_align)):
                if strj_align[idx_i] == '-':
                    curchar   = '-'
                    hashvalue = triplex_score(curchar, stri_align[idx_i], Para)
                else:
                    curchar   = seqtmp[starj + j_idx - 1]
                    hashvalue = triplex_score(curchar, stri_align[idx_i], Para)
                    j_idx += 1
                if curchar == prechar == 'A':
                    tri_score -= prescore
                    tri_score += penaltyT
                    hashvalue  = penaltyT
                if curchar == prechar == 'G':
                    tri_score -= prescore
                    tri_score += penaltyC
                    hashvalue  = penaltyC
                prescore   = hashvalue
                prechar    = curchar
                tri_score += hashvalue

            score      //= 10
            final_score = score / nt
            tri_score  /= nt
            atriplex = Triplex(
                stari=stari, endi=endi,
                starj=starj  + dnaStartPos,
                endj =endj   + dnaStartPos,
                reverse=strand, strand=Para, rule=rule, nt=nt,
                score=final_score, identity=identity, tri_score=tri_score,
                stri_align=stri_align, strj_align=strj_align)

        elif strand == 1 and ntMin <= nt <= ntMax:
            seqtmp = strSrc
            if Para < 0:
                seqtmp = complement(seqtmp)
            j_idx = 0
            hashvalue = prescore = 0.0
            prechar = curchar = ''
            for idx_i in range(len(strj_align)):
                if strj_align[idx_i] == '-':
                    curchar   = '-'
                    hashvalue = triplex_score(curchar, stri_align[idx_i], Para)
                else:
                    curchar   = seqtmp[N - starj - j_idx]
                    hashvalue = triplex_score(curchar, stri_align[idx_i], Para)
                    j_idx += 1
                if curchar == prechar == 'A':
                    tri_score -= prescore + 1000
                    hashvalue  = -1000
                if curchar == prechar == 'G':
                    tri_score -= prescore
                    hashvalue  = 0.0
                prescore   = hashvalue
                prechar    = curchar
                tri_score += hashvalue

            score      //= 10
            final_score = score / nt
            tri_score  /= nt
            atriplex = Triplex(
                stari=stari, endi=endi,
                starj=N - starj + dnaStartPos + 1,
                endj =N - endj  + dnaStartPos + 1,
                reverse=strand, strand=Para, rule=rule, nt=nt,
                score=final_score, identity=identity, tri_score=tri_score,
                stri_align=stri_align, strj_align=strj_align)

        if nt >= ntMin:
            triplex_list.append(atriplex)

        if count == 0:
            break

        # ─── re-expand the DP in the neighbourhood of the current hit ────
        flag = False
        for j in range(nn, n1 - 1, -1):
            CC[j] = 0;  EE[j] = j
            DD[j] = -(Q); FF[j] = j
            RR[j] = SS[j] = mm + 1

        for i in range(mm, m1 - 1, -1):
            c = p = 0
            f  = -(Q)
            ci = fi = i;  pi = i + 1
            cj = fj = pj = nn + 1
            va = V[ord(A[i])]

            for j in range(nn, n1 - 1, -1):
                f -= R
                c2 = c - Q - R
                f, fi, fj = order(f, fi, fj, c2, ci, cj)

                c  = CC[j] - Q - R
                ci = RR[j];  cj = EE[j]
                d  = DD[j] - R
                di = SS[j];  dj = FF[j]
                d, di, dj = order(d, di, dj, c, ci, cj)

                c = 0
                if j not in row[i]:
                    c = p + va[ord(B[j])]

                if c <= 0:
                    c = 0;  ci = i;  cj = j
                else:
                    ci = pi;  cj = pj

                c, ci, cj = order(c, ci, cj, d, di, dj)
                c, ci, cj = order(c, ci, cj, f, fi, fj)

                p  = CC[j];  CC[j] = c
                pi = RR[j];  pj = EE[j]
                RR[j] = ci;  EE[j] = cj
                DD[j] = d
                SS[j] = di;  FF[j] = dj

                if c > 0:
                    flag = True

            HH[i] = CC[n1];  II[i] = RR[n1];  JJ[i] = EE[n1]
            WW[i] = f;       XX[i] = fi;       YY[i] = fj

        rl_ref = [endi - stari + 1]
        cl_ref = [endj - starj + 1]

        while True:
            rflag = cflag = True
            while (rflag and m1 > 1) or (cflag and n1 > 1):

                if rflag and m1 > 1:
                    rflag = False
                    m1 -= 1
                    c = p = 0
                    f  = -(Q)
                    ci = fi = m1;  pi = m1 + 1
                    cj = fj = pj = nn + 1
                    va = V[ord(A[m1])]

                    for j in range(nn, n1 - 1, -1):
                        f -= R
                        c2 = c - Q - R
                        f, fi, fj = order(f, fi, fj, c2, ci, cj)

                        c  = CC[j] - Q - R
                        ci = RR[j];  cj = EE[j]
                        d  = DD[j] - R
                        di = SS[j];  dj = FF[j]
                        d, di, dj = order(d, di, dj, c, ci, cj)

                        c = 0
                        if j not in row[m1]:
                            c = p + va[ord(B[j])]

                        if c <= 0:
                            c = 0;  ci = m1;  cj = j
                        else:
                            ci = pi;  cj = pj

                        c, ci, cj = order(c, ci, cj, d, di, dj)
                        c, ci, cj = order(c, ci, cj, f, fi, fj)

                        p  = CC[j];  CC[j] = c
                        pi = RR[j];  pj = EE[j]
                        RR[j] = ci;  EE[j] = cj
                        DD[j] = d
                        SS[j] = di;  FF[j] = dj

                        if c > 0:
                            flag = True
                        if not rflag and (
                                (ci > rl_ref[0] and cj > cl_ref[0]) or
                                (di > rl_ref[0] and dj > cl_ref[0]) or
                                (fi > rl_ref[0] and fj > cl_ref[0])):
                            rflag = True

                    HH[m1] = CC[n1];  II[m1] = RR[n1];  JJ[m1] = EE[n1]
                    WW[m1] = f;       XX[m1] = fi;       YY[m1] = fj
                    if not cflag and (
                            (ci > rl_ref[0] and cj > cl_ref[0]) or
                            (di > rl_ref[0] and dj > cl_ref[0]) or
                            (fi > rl_ref[0] and fj > cl_ref[0])):
                        cflag = True

                if cflag and n1 > 1:
                    cflag = False
                    n1 -= 1
                    c = 0
                    f  = -(Q)
                    cj = fj = n1
                    va = V[ord(B[n1])]
                    p = 0
                    ci = fi = pi = mm + 1
                    pj = n1 + 1

                    for i in range(mm, m1 - 1, -1):
                        f -= R
                        c2 = c - Q - R
                        f, fi, fj = order(f, fi, fj, c2, ci, cj)

                        c  = HH[i] - Q - R
                        ci = II[i];  cj = JJ[i]
                        d  = WW[i] - R
                        di = XX[i];  dj = YY[i]
                        d, di, dj = order(d, di, dj, c, ci, cj)

                        c = 0
                        if n1 not in row[i]:
                            c = p + va[ord(A[i])]

                        if c <= 0:
                            c = 0;  ci = i;  cj = n1
                        else:
                            ci = pi;  cj = pj

                        c, ci, cj = order(c, ci, cj, d, di, dj)
                        c, ci, cj = order(c, ci, cj, f, fi, fj)

                        p  = HH[i];  HH[i] = c
                        pi = II[i];  pj = JJ[i]
                        II[i] = ci;  JJ[i] = cj
                        WW[i] = d
                        XX[i] = di;  YY[i] = dj

                        if c > 0:
                            flag = True
                        if not cflag and (
                                (ci > rl_ref[0] and cj > cl_ref[0]) or
                                (di > rl_ref[0] and dj > cl_ref[0]) or
                                (fi > rl_ref[0] and fj > cl_ref[0])):
                            cflag = True

                    CC[n1] = HH[m1];  RR[n1] = II[m1];  EE[n1] = JJ[m1]
                    DD[n1] = f;       SS[n1] = fi;       FF[n1] = fj
                    if not rflag and (
                            (ci > rl_ref[0] and cj > cl_ref[0]) or
                            (di > rl_ref[0] and dj > cl_ref[0]) or
                            (fi > rl_ref[0] and fj > cl_ref[0])):
                        rflag = True

            if (m1 == 1 and n1 == 1) or no_cross(LIST, numnode[0], m1, mm, n1, nn, rl_ref, cl_ref):
                break
            m1 -= 1;  n1 -= 1

        if flag:
            for j in range(n1 + 1, nn + 1):
                CC[j] = 0;  RR[j] = m1;  EE[j] = j
                DD[j] = -(Q); SS[j] = m1; FF[j] = j

            for i in range(m1 + 1, mm + 1):
                c = 0;  f = -(Q)
                ci = fi = i;  va = V[ord(A[i])]
                p = 0;  pi = i - 1
                cj = fj = pj = n1

                for j in range(n1 + 1, nn + 1):
                    f -= R
                    c2 = c - Q - R
                    f, fi, fj = order(f, fi, fj, c2, ci, cj)

                    c  = CC[j] - Q - R
                    ci = RR[j];  cj = EE[j]
                    d  = DD[j] - R
                    di = SS[j];  dj = FF[j]
                    d, di, dj = order(d, di, dj, c, ci, cj)

                    c = 0
                    if j not in row[i]:
                        c = p + va[ord(B[j])]

                    if c <= 0:
                        c = 0;  ci = i;  cj = j
                    else:
                        ci = pi;  cj = pj

                    c, ci, cj = order(c, ci, cj, d, di, dj)
                    c, ci, cj = order(c, ci, cj, f, fi, fj)

                    p  = CC[j];  CC[j] = c
                    pi = RR[j];  pj = EE[j]
                    RR[j] = ci;  EE[j] = cj
                    DD[j] = d
                    SS[j] = di;  FF[j] = dj

                    if c > 0:
                        addnode(c, ci, cj, i, j, LIST, numnode)


# ─────────────────────────────────────────────────────────────
# cluster_triplex
# ─────────────────────────────────────────────────────────────

def cluster_triplex(dd: int, length: int,
                    triplex_list: List[Triplex],
                    class1:  List[Dict[int, int]],
                    class1a: List[Dict[int, int]],
                    class1b: List[Dict[int, int]],
                    class_level: int) -> None:
    """
    Cluster triplexes by proximity along the genome axis.
    Fills class1 / class1a / class1b coverage maps.
    """
    axis_map: Dict[int, Axis] = {}
    max_neartriplexnum = 0
    max_pos   = 0
    find      = False
    count     = 0

    for t in triplex_list:
        if t.nt > length:
            count += 1
            middle   = (t.stari + t.endi) // 2
            t.middle = middle
            t.motif  = 0
            if middle not in axis_map:
                axis_map[middle] = Axis()
            axis_map[middle].triplexnum += 1

            for i in range(-dd, dd + 1):
                pos = middle + i
                if pos not in axis_map:
                    axis_map[pos] = Axis()
                if i > 0:
                    axis_map[pos].neartriplex += dd - i
                elif i < 0:
                    axis_map[pos].neartriplex += dd + i

                if axis_map[middle].triplexnum > 0:
                    if axis_map[pos].neartriplex > max_neartriplexnum:
                        max_neartriplexnum = axis_map[pos].neartriplex
                        max_pos = pos
                        find    = True

            t.neartriplex = axis_map[middle].neartriplex

    theclass = 1
    while find:
        for i in range(max_pos - dd, max_pos + dd + 1):
            for t in triplex_list:
                if t.middle == i and t.motif == 0:
                    t.motif  = theclass
                    t.center = max_pos
                    if theclass <= class_level:
                        if t.endj > t.starj:
                            for j in range(t.starj, t.endj):
                                class1 [theclass][j] = class1 [theclass].get(j, 0) + 1
                                class1a[theclass][j] = class1a[theclass].get(j, 0) + 1
                        else:
                            for j in range(t.endj, t.starj):
                                class1 [theclass][j] = class1 [theclass].get(j, 0) + 1
                                class1b[theclass][j] = class1b[theclass].get(j, 0) - 1
            if i in axis_map:
                del axis_map[i]

        max_neartriplexnum = 0
        find = False
        for pos, ax in axis_map.items():
            if ax.neartriplex >= max_neartriplexnum and ax.triplexnum > 0:
                max_neartriplexnum = ax.neartriplex
                max_pos = pos
                find    = True

        theclass += 1


# ─────────────────────────────────────────────────────────────
# print_cluster
# ─────────────────────────────────────────────────────────────

def print_cluster(c_level: int,
                  class1: List[Dict[int, int]],
                  start_genome: int,
                  chro_info: str,
                  dna_size: int,
                  rna_name: str,
                  distance: int,
                  length: int,
                  outFilePath: str,
                  c_tmp_dd: str,
                  c_tmp_length: str,
                  w_tmp_class: List[TmpClass]) -> None:
    """
    Write a UCSC bedGraph track for the cluster at level c_level.
    """
    c_tmp_level = str(c_level)
    class_name  = (outFilePath[:-10] +
                   f"-TFOclass{c_tmp_level}-{c_tmp_dd}-{c_tmp_length}")

    with open(class_name, 'w') as outfile:
        outfile.write(f"browser position {chro_info}:{start_genome}-"
                      f"{start_genome + dna_size}\n")
        outfile.write("browser hide all\n")
        outfile.write("browser pack refGene encodeRegions\n")
        outfile.write("browser full altGraph\n")
        outfile.write("# 300 base wide bar graph, ausoScale is on by default == graphing\n")
        outfile.write("# limits will dynamically change to always show full range of data\n")
        outfile.write("# in viewing window, priority = 20 position this as the second graph\n")
        outfile.write("# Note, zero-relative, half-open coordinate system in use for bedGraph format\n")
        outfile.write(f"track type=bedGraph name='{rna_name} TTS ({c_level})' "
                      f"description='{distance}-{length}' visibility=full "
                      f"color=200,100,0 altColor=0,100,200 priority=20\n")

        data = class1[c_level]
        if not data:
            return

        final_genome = max(data.keys()) + start_genome
        items        = sorted(data.items())

        map_count = map_count1 = 0
        map_tmp1 = map_tmp2 = map_tmp3 = 0
        map_first0 = map_second0 = 0

        idx = 0
        while idx < len(items):
            key0, val0 = items[idx]
            map_tmp1   = key0
            map_tmp2   = val0
            map_first0 = key0
            map_second0 = val0

            if map_count == 0 and map_count1 == 0:
                map_tmp3    = key0
                map_count1 += 1

            if key0 + start_genome == final_genome:
                break

            if idx + 1 >= len(items):
                break
            key1, val1 = items[idx + 1]

            if abs(key1 - map_tmp1) == 1 and val1 == map_tmp2:
                map_tmp1 = key1;  map_tmp2 = val1
                idx += 1
                while idx + 1 < len(items):
                    k2, v2 = items[idx + 1]
                    if abs(k2 - map_tmp1) == 1 and v2 == map_tmp2:
                        map_tmp1 = k2;  map_tmp2 = v2;  idx += 1
                    else:
                        break
                if map_count == 0:
                    w_tmp_class.append(TmpClass(map_first0 + start_genome - 2,
                                                map_tmp1  + start_genome,
                                                map_tmp2, 0, 0))
                    map_count += 1
                else:
                    w_tmp_class.append(TmpClass(map_first0 + start_genome - 1,
                                                map_tmp1  + start_genome,
                                                map_tmp2, 0, 0))
            elif abs(key1 - map_tmp1) != 1 and val1 == map_tmp2:
                w_tmp_class.append(TmpClass(map_tmp1 + start_genome,
                                            key1     + start_genome - 1, 0, 0, 0))
                w_tmp_class.append(TmpClass(key1 + start_genome - 1,
                                            key1 + start_genome, val1, 0, 0))
            elif abs(key1 - map_tmp1) == 1 and val1 != map_tmp2:
                w_tmp_class.append(TmpClass(map_tmp1 + start_genome,
                                            key1     + start_genome, val1, 0, 0))
            elif abs(key1 - map_tmp1) != 1 and val1 != map_tmp2:
                w_tmp_class.append(TmpClass(map_tmp1 + start_genome,
                                            key1     + start_genome - 1, 0, 0, 0))
                w_tmp_class.append(TmpClass(key1 + start_genome - 1,
                                            key1 + start_genome, val1, 0, 0))

            idx += 1

        w_idx = 0
        while w_idx < len(w_tmp_class):
            btc = w_tmp_class[w_idx]
            if btc.genome_start == final_genome:
                break
            if w_idx + 1 == len(w_tmp_class):
                w_idx += 1
                continue
            ctc = w_tmp_class[w_idx + 1]
            if btc.genome_start == ctc.genome_start:
                outfile.write(f"{chro_info}\t{btc.genome_start}\t"
                              f"{ctc.genome_end}\t{ctc.signal_level}\n")
                w_idx += 2
            else:
                outfile.write(f"{chro_info}\t{btc.genome_start}\t"
                              f"{btc.genome_end}\t{btc.signal_level}\n")
                w_idx += 1
