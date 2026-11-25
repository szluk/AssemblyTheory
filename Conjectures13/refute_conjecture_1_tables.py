#!/usr/bin/env python3
"""
verify_criterion.py

Companion to verify_sojas_plain.c for

    S. Lukaszyk, P. Masierak, T. Newton, S. Tworz,
    "A Note on 'Assembly Theory: Optimal Joint Assembly Spaces and Assembly
     Space Energy'", 2026,

refuting Conjecture 1 of

    S. Lukaszyk, P. Masierak, Mathematics 14 (2026), no. 16, 2950,
    doi:10.3390/math14162950                                          [1]

Where verify_sojas_plain.c runs the exhaustive search, this script documents
the counterexample itself and the three statements that surround it in the
Note but are not part of the search:

  1. the incompatibility table of Theorem 1, printed with the intersections
     and both parent pairs of every violating shared term;

  2. the n-ary form of the criterion, validated against the ensembles whose
     SOJAS status is already settled in [1, Proposition 4] -- a regression
     test that the criterion is the right one and not fitted to one pair;

  3. the multiplicity optimal size |mu_T| of the counterexample, and the
     behaviour of the obstruction under scaling.

Criterion (pairs).  For term sets P, Q of two shortest addition chains,
compatible decomposition maps exist if and only if

    for all s in (P cap Q) \\ {1},  exists a, b in P cap Q : s = a + b.

Criterion (n targets).  For shortest chains C_1, ..., C_m, a compatible
assignment exists if and only if every s /= 1 belonging to at least two of the
C_k has a decomposition s = a + b with a, b in the intersection of those C_k
that contain s.

Run
---
    python3 verify_criterion.py

Exit status is 0 if every check passes, 1 otherwise.

SPDX-License-Identifier: MIT
"""

import itertools
import sys

sys.setrecursionlimit(100000)

FAILURES = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'OK ' if ok else 'FAIL'}] {label}: {got}" + ("" if ok else f"  (expected {want})"))
    if not ok:
        FAILURES.append(label)
    return ok


# --------------------------------------------------------------------------
# shortest addition chains: iterative deepening with the doubling bound of
# [1, Appendix A], no further shortcuts
# --------------------------------------------------------------------------

def chains_of_length(N, r):
    out = []

    def rec(c):
        j = len(c) - 1
        if c[-1] == N:
            if j == r:
                out.append(tuple(c))
            return
        if j >= r or c[-1] * (2 ** (r - j)) < N:
            return
        cand = set()
        for a in range(len(c)):
            for b in range(a, len(c)):
                s = c[a] + c[b]
                if c[-1] < s <= N:
                    cand.add(s)
        for s in sorted(cand, reverse=True):
            if s * (2 ** (r - j - 1)) < N:
                continue
            rec(c + [s])

    rec([1])
    return out


_CACHE = {}


def shortest_chains(N):
    """(ASI(N), all shortest addition chains for N)"""
    if N not in _CACHE:
        r = max(1, (N - 1).bit_length())
        while True:
            found = chains_of_length(N, r)
            if found:
                _CACHE[N] = (r, found)
                break
            r += 1
    return _CACHE[N]


def decompositions(s, terms):
    """unordered parent pairs for s inside the term set"""
    t = set(terms)
    return sorted({(a, s - a) for a in t if a <= s - a and (s - a) in t})


def violating_terms(P, Q):
    """shared terms admitting no decomposition inside P cap Q"""
    inter = sorted(set(P) & set(Q))
    return [s for s in inter if s != 1 and not decompositions(s, inter)]


def admits_sojas_pair(M, N):
    for P in shortest_chains(M)[1]:
        for Q in shortest_chains(N)[1]:
            if not violating_terms(P, Q):
                return True
    return False


def admits_sojas(targets):
    """n-ary criterion"""
    families = [shortest_chains(t)[1] for t in targets]
    for combo in itertools.product(*families):
        sets = [set(c) for c in combo]
        ok = True
        for s in sorted(set().union(*sets)):
            if s == 1:
                continue
            owners = [k for k in range(len(sets)) if s in sets[k]]
            if len(owners) < 2:
                continue
            inter = sorted(set.intersection(*[sets[k] for k in owners]))
            if not decompositions(s, inter):
                ok = False
                break
        if ok:
            return True, combo
    return False, None


# --------------------------------------------------------------------------
# 1. the counterexample
# --------------------------------------------------------------------------

def counterexample():
    print("1. The counterexample T = {382, 1387}")

    r382, c382 = shortest_chains(382)
    r1387, c1387 = shortest_chains(1387)
    check("ASI(382)", r382, 11)
    check("ASI(1387)", r1387, 13)
    check("shortest chains for 382", len(c382), 4)
    check("shortest chains for 1387", len(c1387), 4)

    # families differ in the fourth term (8 or not) and in the penultimate term
    fam_p = {c: ("P1" if c[3] == 8 else "P2") for c in c382}
    fam_q = {c: ("Q1" if c[3] == 8 else "Q2") for c in c1387}

    print("\n  Table 1 (all shortest unary chains)")
    for name in ("P1", "P2"):
        members = sorted(c for c in c382 if fam_p[c] == name)
        head = ", ".join(str(v) for v in members[0][:10])
        alts = " or ".join(str(c[10]) for c in sorted(members, key=lambda c: -c[10]))
        print(f"    {name}  {head}, ({alts}), 382")
    for name in ("Q1", "Q2"):
        members = sorted(c for c in c1387 if fam_q[c] == name)
        head = ", ".join(str(v) for v in members[0][:12])
        alts = " or ".join(str(c[12]) for c in sorted(members, key=lambda c: -c[12]))
        print(f"    {name}  {head}, ({alts}), 1387")

    print("\n  incompatibility table (all sixteen pairings)")
    table = {}
    for P in c382:
        for Q in c1387:
            cell = fam_p[P], fam_q[Q]
            inter = sorted(set(P) & set(Q))
            bad = violating_terms(P, Q)
            entry = (tuple(inter), tuple(bad),
                     tuple(decompositions(bad[0], P)) if bad else (),
                     tuple(decompositions(bad[0], Q)) if bad else ())
            table.setdefault(cell, set()).add(entry)

    for cell in sorted(table):
        variants = table[cell]
        if len(variants) != 1:
            FAILURES.append(f"pairing {cell} depends on the penultimate term")
        inter, bad, dp, dq = next(iter(variants))
        pairs_p = " / ".join(str(p) for p in dp)
        pairs_q = " / ".join(str(p) for p in dq)
        print(f"    {cell[0]}/{cell[1]}   P cap Q = {{{', '.join(str(v) for v in inter)}}}"
              f"   violating term {bad[0]}: {pairs_p} in {cell[0]}, {pairs_q} in {cell[1]}")

    check("every pairing has exactly one violating term",
          sorted(next(iter(table[c]))[1][0] for c in sorted(table)), [5, 9, 17, 83])
    check("the pairing outcome is independent of the penultimate term",
          all(len(v) == 1 for v in table.values()), True)
    check("{382, 1387} admits no SOJAS", admits_sojas_pair(382, 1387), False)
    print()


# --------------------------------------------------------------------------
# 2. regression test of the n-ary criterion against [1]
# --------------------------------------------------------------------------

def regression():
    print("2. The n-ary criterion against the ensembles settled in [1]")
    cases = [
        ([17, 23], True, "[1, Section 4]: clash at minimum APU size, SOJAS of size nine"),
        ([43, 77], True, "[1, Proposition 4]: the pair alone admits a SOJAS"),
        ([43, 77, 149], False, "[1, Proposition 4]: smallest unary triple with no SOJAS"),
        ([43, 77, 233], False, "[1, Section 4]"),
        ([43, 77, 382], False, "[1, Section 4]"),
        ([43, 77, 234], True, "control"),
        ([3, 5, 8], True, "[1, Figure 8c]: |pi_T| > |omega_T| but a SOJAS exists"),
    ]
    for targets, expected, note in cases:
        got = admits_sojas_pair(*targets) if len(targets) == 2 else admits_sojas(targets)[0]
        check(f"SOJAS for {targets}  ({note})", got, expected)
    print()


# --------------------------------------------------------------------------
# 3. what survives: MOJAS size, and scaling
# --------------------------------------------------------------------------

def survives():
    print("3. |mu_T| for the counterexample, and the behaviour under scaling")

    c382 = shortest_chains(382)[1]
    c1387 = shortest_chains(1387)[1]

    # Two step vertices are identified in a MOJAS when they carry the same step
    # string AND the same unordered parent pair, [1, Definition 12].  Both plans
    # can realise a shared term s with the same parents exactly when s has a
    # decomposition inside P cap Q, so the maximum overlap is the number of
    # non-violating shared terms.
    best, arg = -1, None
    for P in c382:
        for Q in c1387:
            inter = sorted(set(P) & set(Q))
            good = sum(1 for s in inter if s != 1 and decompositions(s, inter))
            if good > best:
                best, arg = good, (P, Q)
    check("maximum number of identifiable step vertices", best, 3)
    check("|mu_T| = 11 + 13 - 3", 11 + 13 - best, 21)
    print(f"    attained on P cap Q = {{{', '.join(str(v) for v in sorted(set(arg[0]) & set(arg[1])))}}}"
          f", identifying 2, 4, 8 (17 clashes)")

    # |omega_T| >= 14, because no shortest chain for 1387 contains 382
    check("no shortest chain for 1387 contains 382",
          any(382 in c for c in c1387), False)
    print("    hence |omega_T| >= 14, while |pi_T| is undefined")

    print("\n  scaling: the obstruction does not propagate to multiples")
    for mult in (2, 3, 4):
        M, N = 382 * mult, 1387 * mult
        rM, cM = shortest_chains(M)
        rN, cN = shortest_chains(N)
        ok = admits_sojas_pair(M, N)
        print(f"    {{{M}, {N}}}: ASI = {rM}, {rN}; shortest chains {len(cM)}, {len(cN)}; "
              f"SOJAS {'exists' if ok else 'DOES NOT EXIST'}")
        if not ok:
            FAILURES.append(f"unexpected counterexample at {mult}x")
    print()


# --------------------------------------------------------------------------

if __name__ == "__main__":
    print("verify_criterion.py")
    print("=" * 72)
    counterexample()
    print("=" * 72)
    regression()
    print("=" * 72)
    survives()
    print("=" * 72)
    if FAILURES:
        print(f"FAILED: {len(FAILURES)} check(s): " + "; ".join(FAILURES))
        sys.exit(1)
    print("All checks passed.")
