/*
 * verify_sojas_plain.c
 *
 * Independent verification of Theorem 1 of
 *   S. Lukaszyk, P. Masierak, T. Newton, S. Tworz,
 *   "A Note on 'Assembly Theory: Optimal Joint Assembly Spaces and Assembly
 *    Space Energy'", 2026,
 * which refutes Conjecture 1 of
 *   S. Lukaszyk, P. Masierak, Mathematics 14 (2026), no. 16, 2950,
 *   doi:10.3390/math14162950.
 *
 * This is a deliberately plain reimplementation: it uses ONLY the iterative
 * deepening and the doubling bound of [1, Appendix A], and NONE of the
 * last-step / penultimate-step shortcuts of verify_sojas_revised.cpp.  Its
 * purpose is to cross-check the chain counts and the pair test by a route
 * that shares no optimisation with the primary implementation.
 *
 * What it does
 * ------------
 *   1. For every 2 <= N <= B it enumerates ALL shortest addition chains
 *      1 = c_0 < c_1 < ... < c_k = N, starting at k = ceil(log2 N) and
 *      increasing k until a chain is found.  A partial chain ending at c_j is
 *      pruned unless c_j * 2^(k-j) >= N, since no assembly step more than
 *      doubles the current maximum.
 *
 *   2. For every pair of targets 2 <= M < N <= B it tests the compatibility
 *      criterion of Theorem 1,
 *
 *          for all s in (P cap Q) \ {1},  exists a, b in P cap Q : s = a + b,
 *
 *      over all pairs of shortest chains (P for M, Q for N), stopping at the
 *      first compatible pairing.  A pair for which no pairing is compatible
 *      admits no SOJAS.
 *
 *      The criterion constrains only the term sets, not the decomposition
 *      maps, so enumerating chains rather than plans is sufficient.
 *
 * Build
 * -----
 *   gcc -O2 -o verify_sojas_plain verify_sojas_plain.c
 *
 * Run
 * ---
 *   ./verify_sojas_plain [B]            (default B = 1500)
 *
 * Expected output for B = 1500:
 *   3,140,207 shortest chains for 2 <= N <= 1500
 *   1,122,751 pairs tested
 *   exactly one incompatible pair: {382, 1387}
 *
 * Memory: about 100 MB at B = 1500.  Runtime: tens of minutes, dominated by
 * the failed iterative-deepening levels that the shortcuts would have skipped.
 *
 * SPDX-License-Identifier: MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXB 4096       /* hard cap on B                                     */
#define MAXK 24         /* hard cap on chain length                          */

static int B = 1500;

/* ------------------------------------------------------------------ */
/* enumeration of all shortest addition chains                        */
/* ------------------------------------------------------------------ */

static int target, K;
static int cur[MAXK + 1];

static unsigned short *store[MAXB + 1];   /* chains, packed (klen+1) each    */
static int cnt[MAXB + 1];                 /* number of shortest chains       */
static int cap[MAXB + 1];
static int klen[MAXB + 1];                /* ASI(N) = length of the chain    */

static long long total_chains = 0;

static void push_chain(void)
{
    int N = target, k = K;

    if (cnt[N] == cap[N]) {
        cap[N] = cap[N] ? cap[N] * 2 : 8;
        store[N] = (unsigned short *)realloc(
                       store[N], (size_t)cap[N] * (k + 1) * sizeof(unsigned short));
        if (!store[N]) { fprintf(stderr, "out of memory at N = %d\n", N); exit(1); }
    }
    unsigned short *p = store[N] + (size_t)cnt[N] * (k + 1);
    for (int j = 0; j <= k; j++) p[j] = (unsigned short)cur[j];
    cnt[N]++;
}

static char seen_sum[2 * MAXB + 4];

/* cur[0..j] is a strictly increasing partial chain; extend it. */
static void rec(int j)
{
    if (cur[j] == target) { if (j == K) push_chain(); return; }
    if (j >= K) return;

    int rem = K - j;                       /* steps still available          */

    /* doubling bound on the partial chain itself */
    long long mx = cur[j];
    for (int t = 0; t < rem; t++) { mx <<= 1; if (mx >= target) break; }
    if (mx < target) return;

    /* distinct sums of two earlier terms, strictly above the current maximum */
    int sums[(MAXK + 1) * (MAXK + 2) / 2], ns = 0;
    for (int a = 0; a <= j; a++)
        for (int b = a; b <= j; b++) {
            int s = cur[a] + cur[b];
            if (s > cur[j] && s <= target && !seen_sum[s]) {
                seen_sum[s] = 1;
                /* doubling bound on the candidate */
                long long v = s;
                for (int t = 0; t < rem - 1; t++) { v <<= 1; if (v >= target) break; }
                if (v >= target) sums[ns++] = s;
            }
        }
    for (int a = 0; a <= j; a++)           /* clear the markers              */
        for (int b = a; b <= j; b++) {
            int s = cur[a] + cur[b];
            if (s > cur[j] && s <= target) seen_sum[s] = 0;
        }

    for (int t = ns - 1; t >= 0; t--) { cur[j + 1] = sums[t]; rec(j + 1); }
}

static int ceil_log2(int n) { int k = 0; while ((1 << k) < n) k++; return k; }

static void enumerate(int N)
{
    target = N;
    if (N == 1) { K = 0; cur[0] = 1; push_chain(); klen[1] = 0; return; }
    for (K = ceil_log2(N); K <= MAXK; K++) {
        cur[0] = 1;
        rec(0);
        if (cnt[N]) { klen[N] = K; return; }
    }
    fprintf(stderr, "no chain found for %d\n", N);
    exit(1);
}

/* ------------------------------------------------------------------ */
/* compatibility criterion                                            */
/* ------------------------------------------------------------------ */

static char inP[MAXB + 2];                 /* membership mask of P           */

/* Every shared term other than 1 must be the sum of two shared terms. */
static int compatible(const unsigned short *Q, int kq)
{
    int I[MAXK + 2], ni = 0;
    for (int j = 0; j <= kq; j++) if (inP[Q[j]]) I[ni++] = Q[j];   /* sorted  */

    for (int a = 0; a < ni; a++) {
        int s = I[a];
        if (s == 1) continue;
        int ok = 0;
        for (int b = 0; b < ni && !ok; b++) {
            if (I[b] * 2 > s) break;       /* need I[b] <= s/2               */
            for (int c = b; c < ni; c++) {
                if (I[b] + I[c] == s) { ok = 1; break; }
                if (I[b] + I[c] >  s) break;
            }
        }
        if (!ok) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    if (argc > 1) B = atoi(argv[1]);
    if (B < 2 || B > MAXB) { fprintf(stderr, "B must be in [2, %d]\n", MAXB); return 1; }

    for (int N = 1; N <= B; N++) { enumerate(N); total_chains += cnt[N]; }

    printf("B = %d\n", B);
    printf("shortest chains, 1 <= N <= B : %lld\n", total_chains);
    printf("shortest chains, 2 <= N <= B : %lld\n", total_chains - cnt[1]);
    if (B >= 382)  printf("ASI(382)  = %d  (%d shortest chains)\n", klen[382],  cnt[382]);
    if (B >= 1387) printf("ASI(1387) = %d  (%d shortest chains)\n", klen[1387], cnt[1387]);
    fflush(stdout);

    long long npairs = 0, nbad = 0;
    for (int M = 2; M <= B; M++) {
        for (int N = M + 1; N <= B; N++) {
            npairs++;
            int found = 0;
            for (int a = 0; a < cnt[M] && !found; a++) {
                const unsigned short *P = store[M] + (size_t)a * (klen[M] + 1);
                for (int j = 0; j <= klen[M]; j++) inP[P[j]] = 1;
                for (int b = 0; b < cnt[N]; b++) {
                    const unsigned short *Q = store[N] + (size_t)b * (klen[N] + 1);
                    if (compatible(Q, klen[N])) { found = 1; break; }
                }
                for (int j = 0; j <= klen[M]; j++) inP[P[j]] = 0;
            }
            if (!found) {
                nbad++;
                printf("NO SOJAS: {%d, %d}\n", M, N);
                fflush(stdout);
            }
        }
        if (M % 100 == 0) fprintf(stderr, "M = %d\n", M);
    }

    printf("pairs tested: %lld, incompatible: %lld\n", npairs, nbad);
    return 0;
}
