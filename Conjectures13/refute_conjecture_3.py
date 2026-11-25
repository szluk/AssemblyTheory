#!/usr/bin/env python3
"""
verify_esojas.py

Independent verification of Theorem 2 of

    S. Lukaszyk, P. Masierak, T. Newton, S. Tworz,
    "A Note on 'Assembly Theory: Optimal Joint Assembly Spaces and Assembly
     Space Energy'", 2026,

which refutes Conjecture 3 of

    S. Lukaszyk, P. Masierak, Mathematics 14 (2026), no. 16, 2950,
    doi:10.3390/math14162950                                          [1]

The script establishes, by exhaustive enumeration, the two claims that the
finite lists in the Note cannot establish by themselves, namely that Table 2
and Table 3 are COMPLETE lists of the shortest plans, and that the minimum
collective sizes are 5 and 7.  Everything else in Theorem 2 (the energies, the
union sizes, the clash) follows from those lists by arithmetic, and is
recomputed here as well.

Definitions used
----------------
An assembly plan for a target is a set of step strings, each equal to the
ordered concatenation of two strings that are either basic symbols or other
step strings of the plan.  Parents are strictly shorter than their child, so a
plan is acyclic and the step strings of a minimum-size plan are distinct.

Pathway depth and Dirichlet energy, [1, Eqs. (2)-(4) and Proposition 1]:

    d(s) = 1 + max{d(u), d(v)},   d(basic symbol) = 0,
    E(H) = sum over step vertices s of [ 1 + (1 + |d(u) - d(v)|)^2 ].

Two plans are compatible (their union is an SOJAS candidate) when every step
string occurring in both carries the same parent pair, the pair being taken
unordered, as in [1, Definition 12].

Run
---
    python3 verify_esojas.py

Exit status is 0 if every claim of Theorem 2 is reproduced, 1 otherwise.

SPDX-License-Identifier: MIT
"""

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
# general string case
# --------------------------------------------------------------------------

def splits(s):
    """all ordered decompositions s = a o b"""
    return [(s[:j], s[j:]) for j in range(1, len(s))]


def joint_plans(targets, k):
    """all plans of exactly k steps whose step strings cover every target"""
    out = []

    def rec(assigned, pending):
        if len(assigned) + len(pending) > k:
            return
        if not pending:
            if len(assigned) == k:
                out.append(dict(assigned))
            return
        s = max(pending, key=len)
        rest = pending - {s}
        for a, b in splits(s):
            nxt = dict(assigned)
            nxt[s] = (a, b)
            pend = set(rest)
            for x in (a, b):
                if len(x) > 1 and x not in nxt:
                    pend.add(x)
            rec(nxt, pend)

    rec({}, set(targets))
    seen, uniq = set(), []
    for p in out:
        key = tuple(sorted((s, a, b) for s, (a, b) in p.items()))
        if key not in seen:
            seen.add(key)
            uniq.append(p)
    return uniq


def min_joint(targets, kmax=12):
    """(minimum number of steps, all plans attaining it)"""
    for k in range(1, kmax + 1):
        plans = joint_plans(targets, k)
        if plans:
            return k, plans
    raise RuntimeError("no plan found within kmax")


def depth_map(plan):
    d = {}

    def dd(x):
        if len(x) == 1:
            return 0
        if x not in d:
            a, b = plan[x]
            d[x] = 1 + max(dd(a), dd(b))
        return d[x]

    for s in plan:
        dd(s)
    return d


def energy(plan):
    d = depth_map(plan)
    tot = 0
    for s, (a, b) in plan.items():
        da = 0 if len(a) == 1 else d[a]
        db = 0 if len(b) == 1 else d[b]
        tot += 1 + (1 + abs(da - db)) ** 2
    return tot


def union(p, q):
    """(union, None) if compatible, (None, clashing string) otherwise"""
    u = dict(p)
    for s, ab in q.items():
        if s in u and tuple(sorted(u[s])) != tuple(sorted(ab)):
            return None, s
        u[s] = ab
    return u, None


def show_plan(plan):
    order = sorted(plan, key=lambda s: (len(s), s))
    return ", ".join(f"{s} = {plan[s][0]} o {plan[s][1]}" for s in order)


# --------------------------------------------------------------------------
# unary case (strings identified with their lengths)
# --------------------------------------------------------------------------

def joint_plans_u(targets, k):
    out = []

    def rec(assigned, pending):
        if len(assigned) + len(pending) > k:
            return
        if not pending:
            if len(assigned) == k:
                out.append(dict(assigned))
            return
        s = max(pending)
        rest = pending - {s}
        for a in range(1, s // 2 + 1):          # unordered parent pairs
            b = s - a
            nxt = dict(assigned)
            nxt[s] = (a, b)
            pend = set(rest)
            for x in (a, b):
                if x > 1 and x not in nxt:
                    pend.add(x)
            rec(nxt, pend)

    rec({}, set(targets))
    seen, uniq = set(), []
    for p in out:
        key = tuple(sorted((s, a, b) for s, (a, b) in p.items()))
        if key not in seen:
            seen.add(key)
            uniq.append(p)
    return uniq


def min_joint_u(targets, kmax=12):
    for k in range(1, kmax + 1):
        plans = joint_plans_u(targets, k)
        if plans:
            return k, plans
    raise RuntimeError("no plan found within kmax")


def energy_u(plan):
    d = {}

    def dd(x):
        if x == 1:
            return 0
        if x not in d:
            a, b = plan[x]
            d[x] = 1 + max(dd(a), dd(b))
        return d[x]

    for s in plan:
        dd(s)
    tot = 0
    for s, (a, b) in plan.items():
        da = 0 if a == 1 else d[a]
        db = 0 if b == 1 else d[b]
        tot += 1 + (1 + abs(da - db)) ** 2
    return tot


def union_u(p, q):
    u = dict(p)
    for s, ab in q.items():
        if s in u and tuple(sorted(u[s])) != tuple(sorted(ab)):
            return None, s
        u[s] = ab
    return u, None


def key_u(plan, omit=()):
    return tuple(sorted((s, min(ab), max(ab)) for s, ab in plan.items() if s not in omit))


# --------------------------------------------------------------------------
# Theorem 2, binary case: T = {00001, 00101}
# --------------------------------------------------------------------------

TABLE2 = {
    "A1": (("00", "0", "0"),   ("001",  "00", "1"),  ("00001", "00",   "001")),
    "A2": (("00", "0", "0"),   ("0000", "00", "00"), ("00001", "0000", "1")),
    "B1": (("01", "0", "1"),   ("0101", "01", "01"), ("00101", "0",    "0101")),
    "B2": (("01", "0", "1"),   ("001",  "0",  "01"), ("00101", "001",  "01")),
}


def binary_case():
    x, y = "00001", "00101"
    print(f"Binary case, T = {{{x}, {y}}}")

    kx, px = min_joint([x])
    ky, py = min_joint([y])
    check(f"ASI({x})", kx, 3)
    check(f"ASI({y})", ky, 3)
    check(f"number of shortest plans for {x}", len(px), 2)
    check(f"number of shortest plans for {y}", len(py), 2)

    # match the enumerated plans against Table 2
    def as_key(plan):
        return tuple(sorted((s, a, b) for s, (a, b) in plan.items()))

    label_of = {}
    for name, rows in TABLE2.items():
        label_of[tuple(sorted(rows))] = name

    found = {}
    for plan in px + py:
        k = as_key(plan)
        if k in label_of:
            found[label_of[k]] = plan
    check("Table 2 reproduced exactly", sorted(found), ["A1", "A2", "B1", "B2"])

    print("\n  Table 2 (all shortest binary plans)")
    print(f"    {'Plan':5s} {'steps':48s} {'E':>3s}")
    for name in ("A1", "A2", "B1", "B2"):
        p = found[name]
        print(f"    {name:5s} {show_plan(p):48s} {energy(p):>3d}")
    check("energies (A1, A2, B1, B2)",
          tuple(energy(found[n]) for n in ("A1", "A2", "B1", "B2")), (12, 14, 14, 12))

    print("\n  unions")
    sizes, clashes = {}, {}
    for a in ("A1", "A2"):
        for b in ("B1", "B2"):
            u, cl = union(found[a], found[b])
            if u is None:
                clashes[(a, b)] = cl
                print(f"    {a} u {b}: CLASH at {cl}")
            else:
                sizes[(a, b)] = (len(u), energy(u))
                print(f"    {a} u {b}: size {len(u)}, E = {energy(u)}")
    check("the only clash", sorted(clashes.items()), [(("A1", "B2"), "001")])
    check("union energies (A1uB1, A2uB1, A2uB2)",
          (sizes[("A1", "B1")][1], sizes[("A2", "B1")][1], sizes[("A2", "B2")][1]),
          (26, 28, 26))

    pi = min(v[0] for v in sizes.values())
    sojas = {k: v for k, v in sizes.items() if v[0] == pi}
    emin = min(v[1] for v in sojas.values())
    esojas = [k for k, v in sojas.items() if v[1] == emin]
    kw, pw = min_joint([x, y])

    check("|pi_T|", pi, 6)
    check("|omega_T|", kw, 5)
    check("E(pihat_T)", emin, 26)
    check("number of ESOJAS", len(esojas), 2)

    print("\n  ESOJAS step-string sets (for the non-isomorphism claim)")
    for a, b in esojas:
        u, _ = union(found[a], found[b])
        print(f"    {a} u {b}: {{" + ", ".join(sorted(u, key=lambda s: (len(s), s))) + "}")
    s1, _ = union(found[esojas[0][0]], found[esojas[0][1]])
    s2, _ = union(found[esojas[1][0]], found[esojas[1][1]])
    check("the two ESOJAS carry different step-string sets", set(s1) != set(s2), True)

    print("\n  a five-step collective construction")
    print("    " + show_plan(pw[0]))
    print()


# --------------------------------------------------------------------------
# Theorem 2, unary case: T = {7, 34}
# --------------------------------------------------------------------------

TABLE3_P = [
    ((3, 1, 2), (4, 1, 3), (7, 3, 4)),
    ((3, 1, 2), (4, 2, 2), (7, 3, 4)),
    ((3, 1, 2), (5, 2, 3), (7, 2, 5)),
    ((3, 1, 2), (6, 3, 3), (7, 1, 6)),
    ((4, 2, 2), (5, 1, 4), (7, 2, 5)),
    ((4, 2, 2), (6, 2, 4), (7, 1, 6)),
]
TABLE3_Q = [
    ((4, 2, 2), (8, 4, 4), (9, 1, 8),   (17, 8, 9),   (34, 17, 17)),
    ((4, 2, 2), (8, 4, 4), (16, 8, 8),  (17, 1, 16),  (34, 17, 17)),
    ((4, 2, 2), (8, 4, 4), (16, 8, 8),  (18, 2, 16),  (34, 16, 18)),
    ((4, 2, 2), (8, 4, 4), (16, 8, 8),  (32, 16, 16), (34, 2, 32)),
]


def unary_case():
    m, n = 7, 34
    print(f"Unary case, T = {{{m}, {n}}}")

    km, pm = min_joint_u([m])
    kn, pn = min_joint_u([n])
    check(f"ASI({m})", km, 4)
    check(f"ASI({n})", kn, 6)
    check(f"number of shortest plans for {m}", len(pm), 6)
    check(f"number of shortest plans for {n}", len(pn), 4)

    # Table 3 omits the common step 2 = 1 + 1; sort lexicographically on the
    # remaining steps, which is the order printed in the Note.
    pm = sorted(pm, key=lambda p: key_u(p, omit=(2,)))
    pn = sorted(pn, key=lambda p: key_u(p, omit=(2,)))
    check("Table 3, plans for 7 reproduced exactly",
          [key_u(p, omit=(2,)) for p in pm], [tuple(r) for r in TABLE3_P])
    check("Table 3, plans for 34 reproduced exactly",
          [key_u(p, omit=(2,)) for p in pn], [tuple(r) for r in TABLE3_Q])

    print("\n  Table 3 (all shortest unary plans; the common step 2 = 1 + 1 omitted)")
    for idx, p in enumerate(pm, 1):
        steps = ", ".join(f"{s} = {a} + {b}" for s, (a, b) in sorted(p.items()) if s != 2)
        print(f"    P{idx}  {steps}")
    for idx, q in enumerate(pn, 1):
        steps = ", ".join(f"{s} = {a} + {b}" for s, (a, b) in sorted(q.items()) if s != 2)
        print(f"    Q{idx}  {steps}")

    check("every shortest plan for 34 contains 2 = 1+1 and 4 = 2+2",
          all(q.get(2) == (1, 1) and q.get(4) == (2, 2) for q in pn), True)

    rows, clashes = {}, []
    for a, p in enumerate(pm, 1):
        for b, q in enumerate(pn, 1):
            u, cl = union_u(p, q)
            if u is None:
                clashes.append((f"P{a}", f"Q{b}", cl))
            else:
                rows[(a, b)] = (len(u), energy_u(u))
    check("clashes occur only for P1 (which has 4 = 1 + 3), always at 4",
          sorted({(c[0], c[2]) for c in clashes}), [("P1", 4)])

    pi = min(v[0] for v in rows.values())
    check("|pi_T|", pi, 8)
    eight = sorted(k for k, v in rows.items() if v[0] == 8)
    check("size-eight unions are exactly P2, P5, P6 against Q1..Q4",
          sorted({f"P{a}" for a, _ in eight}), ["P2", "P5", "P6"])

    print("\n  energies of the twelve size-eight unions")
    print("           " + "".join(f"Q{b:<5d}" for b in range(1, 5)))
    for a in (2, 5, 6):
        print(f"    P{a}     " + "".join(f"{rows[(a, b)][1]:<6d}" for b in range(1, 5)))

    sojas = {k: v for k, v in rows.items() if v[0] == pi}
    emin = min(v[1] for v in sojas.values())
    esojas = sorted(k for k, v in sojas.items() if v[1] == emin)
    check("E(pihat_T)", emin, 37)
    check("number of ESOJAS", len(esojas), 2)
    check("the two ESOJAS", [f"P{a} u Q{b}" for a, b in esojas], ["P2 u Q1", "P2 u Q3"])

    print("\n  ESOJAS step-string sets (for the non-isomorphism claim)")
    sets = []
    for a, b in esojas:
        u, _ = union_u(pm[a - 1], pn[b - 1])
        sets.append(set(u))
        print(f"    P{a} u Q{b}: {{" + ", ".join(str(v) for v in sorted(u)) + "}")
    check("the two ESOJAS carry different step-string sets", sets[0] != sets[1], True)

    kw, pw = min_joint_u([m, n])
    check("|omega_T|", kw, 7)
    print("\n  a seven-step collective construction")
    print("    " + ", ".join(f"{s} = {a} + {b}" for s, (a, b) in sorted(pw[0].items())))
    print()


# --------------------------------------------------------------------------

if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print("=" * 72)
    binary_case()
    print("=" * 72)
    unary_case()
    print("=" * 72)
    if FAILURES:
        print(f"FAILED: {len(FAILURES)} check(s): " + "; ".join(FAILURES))
        sys.exit(1)
    print("All claims of Theorem 2 reproduced.")
