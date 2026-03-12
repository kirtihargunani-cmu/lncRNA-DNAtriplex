"""
rules.py – Python translation of rules.txt

Provides:
  • complement(seq)       – in-place-style DNA complement (returns new string)
  • reverse_seq(seq)      – reverse a sequence string
  • transfer_string(...)  – map a nucleotide sequence through a triplex rule
"""

from __future__ import annotations
import sys

# ─────────────────────────────────────────────────────────────────────────────
# Triplex rule patterns
# Each 10-character string encodes:
#   characters 0-4  → the five input bases to match
#   characters 5-9  → the corresponding output bases
# ─────────────────────────────────────────────────────────────────────────────

# Parallel rules (Para >= 0, strand == 0)
PARARULE1    = "ATGCNTGGTN"
PARARULE2    = "ATGCNTGCTN"
PARARULE3    = "ATGCNTGTTN"
PARARULE4    = "ATGCNTGGCN"
PARARULE5    = "ATGCNTGCCN"
PARARULE6    = "ATGCNTGTCN"

# Parallel rules reversed (Para >= 0, strand == 1)
PARARULE1REV = "ATGCNGTTGN"
PARARULE2REV = "ATGCNGTTCN"
PARARULE3REV = "ATGCNGTTAN"  # note: C uses "ATGCNGTTTN" — kept from source
PARARULE3REV = "ATGCNGTTTN"
PARARULE4REV = "ATGCNGTCGN"
PARARULE5REV = "ATGCNGTCCN"
PARARULE6REV = "ATGCNGTCTN"

# Anti-parallel rules (Para < 0, strand == 0)
ANTIRULE1    = "ATGCNGTTGN"
ANTIRULE2    = "ATGCNGTTCN"
ANTIRULE3    = "ATGCNGTTAN"
ANTIRULE4    = "ATGCNGTCGN"
ANTIRULE5    = "ATGCNGTCCN"
ANTIRULE6    = "ATGCNGTCAN"
ANTIRULE7    = "ATGCNGATGN"
ANTIRULE8    = "ATGCNGATCN"
ANTIRULE9    = "ATGCNGATAN"
ANTIRULE10   = "ATGCNGACGN"
ANTIRULE11   = "ATGCNGACCN"
ANTIRULE12   = "ATGCNGACAN"
ANTIRULE13   = "ATGCNGCTGN"
ANTIRULE14   = "ATGCNGCTCN"
ANTIRULE15   = "ATGCNGCTAN"
ANTIRULE16   = "ATGCNGCCGN"
ANTIRULE17   = "ATGCNGCCCN"
ANTIRULE18   = "ATGCNGCCAN"

# Anti-parallel rules reversed (Para < 0, strand == 1)
ANTIRULE1REV  = "ATGCNTGGTN"
ANTIRULE2REV  = "ATGCNTGCTN"
ANTIRULE3REV  = "ATGCNTGATN"
ANTIRULE4REV  = "ATGCNTGGCN"
ANTIRULE5REV  = "ATGCNTGCCN"
ANTIRULE6REV  = "ATGCNTGACN"
ANTIRULE7REV  = "ATGCNAGGTN"
ANTIRULE8REV  = "ATGCNAGCTN"
ANTIRULE9REV  = "ATGCNAGATN"
ANTIRULE10REV = "ATGCNAGGCN"
ANTIRULE11REV = "ATGCNAGCCN"
ANTIRULE12REV = "ATGCNAGACN"
ANTIRULE13REV = "ATGCNCGGTN"
ANTIRULE14REV = "ATGCNCGCTN"
ANTIRULE15REV = "ATGCNCGATN"
ANTIRULE16REV = "ATGCNCGGCN"
ANTIRULE17REV = "ATGCNCGCCN"
ANTIRULE18REV = "ATGCNCGACN"

# ── Lookup tables for quick rule selection ────────────────────────────────────

_PARA_FWD: dict[int, str] = {
    1: PARARULE1, 2: PARARULE2, 3: PARARULE3,
    4: PARARULE4, 5: PARARULE5, 6: PARARULE6,
}

_PARA_REV: dict[int, str] = {
    1: PARARULE1REV, 2: PARARULE2REV, 3: PARARULE3REV,
    4: PARARULE4REV, 5: PARARULE5REV, 6: PARARULE6REV,
}

_ANTI_FWD: dict[int, str] = {
    1:  ANTIRULE1,  2:  ANTIRULE2,  3:  ANTIRULE3,
    4:  ANTIRULE4,  5:  ANTIRULE5,  6:  ANTIRULE6,
    7:  ANTIRULE7,  8:  ANTIRULE8,  9:  ANTIRULE9,
    10: ANTIRULE10, 11: ANTIRULE11, 12: ANTIRULE12,
    13: ANTIRULE13, 14: ANTIRULE14, 15: ANTIRULE15,
    16: ANTIRULE16, 17: ANTIRULE17, 18: ANTIRULE18,
}

_ANTI_REV: dict[int, str] = {
    1:  ANTIRULE1REV,  2:  ANTIRULE2REV,  3:  ANTIRULE3REV,
    4:  ANTIRULE4REV,  5:  ANTIRULE5REV,  6:  ANTIRULE6REV,
    7:  ANTIRULE7REV,  8:  ANTIRULE8REV,  9:  ANTIRULE9REV,
    10: ANTIRULE10REV, 11: ANTIRULE11REV, 12: ANTIRULE12REV,
    13: ANTIRULE13REV, 14: ANTIRULE14REV, 15: ANTIRULE15REV,
    16: ANTIRULE16REV, 17: ANTIRULE17REV, 18: ANTIRULE18REV,
}


# ─────────────────────────────────────────────────────────────────────────────
# complement
# ─────────────────────────────────────────────────────────────────────────────

_COMP_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")

def complement(seq: str) -> str:
    """
    Return the Watson-Crick complement of a DNA sequence.
    'N' maps to 'N'; any other character is dropped (matches C behaviour).
    """
    result: list[str] = []
    for ch in seq:
        if   ch == 'A': result.append('T')
        elif ch == 'C': result.append('G')
        elif ch == 'G': result.append('C')
        elif ch == 'T': result.append('A')
        elif ch == 'N': result.append('N')
        # unrecognised characters are silently skipped (matches C default:break)
    return ''.join(result)


# ─────────────────────────────────────────────────────────────────────────────
# reverse_seq
# ─────────────────────────────────────────────────────────────────────────────

def reverse_seq(seq: str) -> str:
    """Return the reverse of a sequence string."""
    return seq[::-1]


# ─────────────────────────────────────────────────────────────────────────────
# transfer_string
# ─────────────────────────────────────────────────────────────────────────────

def transfer_string(seq1: str, strand: int, Para: int, rule: int) -> str:
    """
    Map each base in seq1 to a new base using the appropriate triplex rule.

    Parameters
    ----------
    seq1   : input nucleotide sequence
    strand : 0 = forward strand, 1 = reverse strand
    Para   : >= 0 → parallel rules;  < 0 → anti-parallel rules
    rule   : rule number (1-6 for parallel, 1-18 for anti-parallel)

    Returns
    -------
    Transformed sequence of the same length as seq1.
    Unknown/unmatched bases are replaced with 'N'.

    Raises
    ------
    SystemExit (matching C exit(1)) if the rule is not found or the output
    length does not match the input length.
    """
    # ── select rule string ────────────────────────────────────────────────
    if Para >= 0:
        table = _PARA_FWD if strand == 0 else _PARA_REV
    else:
        table = _ANTI_FWD if strand == 0 else _ANTI_REV

    rule_seq = table.get(rule)
    if rule_seq is None:
        sys.exit(1)

    # rule_seq layout:
    #   indices 0-4  → input  bases (keys)
    #   indices 5-9  → output bases (values)
    in_bases  = rule_seq[0:5]   # e.g. "ATGCN"
    out_bases = rule_seq[5:10]  # e.g. "TGGTN"

    # build a per-character mapping dict for O(1) lookups
    mapping: dict[str, str] = {in_bases[k]: out_bases[k] for k in range(5)}

    # ── apply the mapping ─────────────────────────────────────────────────
    tmp_seq = ''.join(mapping.get(ch, 'N') for ch in seq1)

    if len(tmp_seq) != len(seq1):   # should never happen, mirrors C guard
        sys.exit(1)

    return tmp_seq
