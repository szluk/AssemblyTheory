# Assembly Theory: Formalizing Assembly Spaces, Discovering Patterns and Bounds
> Wawrzyniec Bieniawski, Szymon Łukaszyk, Piotr Masierak, Andrzej Tomski, Szymon Tworz

> https://doi.org/10.1098/rsos.260082 (May 2026)

# Assembly Theory: Optimal Joint Assembly Spaces and Assembly Space Energy
> Szymon Łukaszyk, Piotr Masierak

> https://doi.org/10.20944/preprints202607.2172.v1 (June 2026)

> b - alphabet size, N - string length, ASI - assembly index,

> Version: v1 (30.03.2026) 

> Version: v2 (25.06.2026) - including the second paper

## List of files:
---
* A003313.xlsx
> Length of shortest addition chain for N (OEIS A003313) and the minimum assembly index of strings for all b, illustrating particular pathways shown in Tables 1 & 2 and discussed in Section 3 of the paper.
---
* Plans_b2
> List of assembly plans of max ASI strings for b=2 (Tables 3 & 4 of the paper).
---
* Plans_b3
> List of assembly plans of max ASI strings for b=3  (Table 5 of the paper).
---
* Plans_b4
> List of assembly plans of max ASI strings for b=4  (Table 6 of the paper).
---
* Figure10.xls
> Data for Figure 10 and Table 7 of the paper:
"The ASI (a), Constrained ASI (b), Re-Pair steps (c), and LZW steps (d) of all bitstrings of length N=10, as a function of their decimal representations ({0, 1, . . . , 1023}), in squares with sides equal to 32 pixels..." 
---
* allN_2-16.xls
> Lists of all bitstrings of length 2<=N<=16 along with their ASI.
---
* QuickASISolver.exe
> A Windows console program to calculate the string ASI.
---
* SASI-MASI-CASI_N8_b2.xlsx
> Lists of all ensembles of two bitstrings of length N=8 along with their ASI, SASI, MASI, CASI, Dirichlet energy and other relevant parameters.
(SASI-CASI_N8_b2.xlsx is a previous version of this file)
---
* SASI-MASI-CASI_N2-30_b1.md
> Complete enumeration of all COJAS ω_T and all SOJAS π_T for every pair of unary strings T = {m, n} with |ω_T| < |π_T| for 2 ≤ m < n ≤ 30, including the Dirichlet energy of every space.
(SASI-CASI_N2-30_b1.md is a previous version of this file)
---
* minimum_ASI_strings_by_plan.md
> Distinct minimum-ASI strings of length N over an alphabet containing b symbols, broken down by individual shortest addition chain (plan). The last row of each block (in bold) is the union over all plans of that N. The final column gives the count for b = 3.
---
* delta_canonical_N2-148.md
> A global decomposition table: for every unary target 2 ≤ N ≤ 148 it gives one unordered decomposition N = a + b. Setting cl(1) = {} and cl(N) = {N} ∪ cl(a) ∪ cl(b), every closure satisfies |cl(N)| = ASI(N) and is therefore an optimal assembly plan of N; this is verifiable from the table alone, in one pass. Since all plans are read off a single table, no step string ever receives two different decompositions, so the closures of the targets of any ensemble whose targets are all at most 148 form a conflict-free choice of individually optimal plans, whose union is a SOJAS.
