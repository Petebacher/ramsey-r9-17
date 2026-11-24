/* Exact checker by plain enumeration (2026-11-24). Written from scratch; shares no code with check_russian_doll.c, check_split_branch_bound.c or check_branch_bound.c.
   Method: plain exhaustive enumeration, with NO pruning bound of any kind.
   It visits every clique (mode c) or every independent set (mode i) of size < T, in increasing vertex order,
   and counts the sets of every size 1..T. "No T-set" is proved by count[T] == 0.
   The counts of the smaller sizes are a fingerprint that can be cross-checked (count[2] = edges, count[3] = triangles).
   Usage: check_enumerate file T c|i     (file: 0/1 matrix, n rows; spaces allowed; every line must be a row of n entries)
   Parallel over the first vertex (OpenMP). Memory: a few MB. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
#define MAXN 1600
#define MW (MAXN/64)
static int n, nw, T;
static uint64_t g[MAXN][MW];           /* the graph searched for cliques (complement in mode i) */

static void rec(const uint64_t *C, int lo, int size, unsigned long long *cnt) {
  /* C = common neighbours of the current set (size vertices) that have a larger index than all of them.
     words below lo are not valid and not read. */
  uint64_t D[MW], C2[MW];
  for (int i = lo; i < nw; i++) D[i] = C[i];
  for (int w = lo; w < nw; w++) {
    while (D[w]) {
      int b = __builtin_ctzll(D[w]);
      D[w] &= D[w] - 1;                 /* remove v; D now holds only candidates larger than v */
      int v = w * 64 + b;
      unsigned long long p = 0;
      for (int i = w; i < nw; i++) { C2[i] = D[i] & g[v][i]; p += __builtin_popcountll(C2[i]); }
      cnt[size + 2] += p;               /* each vertex of C2 extends (set + v) to a set of size+2 */
      if (p && size + 2 < T) rec(C2, w, size + 1, cnt);
    }
  }
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: check_enumerate file T c|i\n"); return 2; }
  T = atoi(argv[2]); int ind = argv[3][0] == 'i';
  if (argv[3][0] != 'i' && argv[3][0] != 'c') { fprintf(stderr, "mode must be c or i\n"); return 2; }
  if (T < 2 || T > 30) { fprintf(stderr, "T out of range\n"); return 2; }
  FILE *f = fopen(argv[1], "r"); if (!f) { perror("open"); return 2; }
  /* strict parse: count entries per line, all lines must have the same count n, and there must be n lines */
  static unsigned char a[MAXN][MAXN];
  int row = 0, col = 0, ch; n = -1;
  while ((ch = fgetc(f)) != EOF) {
    if (ch == '0' || ch == '1') {
      if (row >= MAXN || col >= MAXN) { fprintf(stderr, "ERROR: larger than MAXN=%d\n", MAXN); return 2; }
      a[row][col++] = ch - '0';
    } else if (ch == '\n') {
      if (col == 0) continue;
      if (n < 0) n = col;
      if (col != n) { fprintf(stderr, "ERROR: row %d has %d entries, expected %d\n", row, col, n); return 2; }
      row++; col = 0;
    } else if (ch != ' ' && ch != '\r') { fprintf(stderr, "ERROR: unexpected character %d\n", ch); return 2; }
  }
  fclose(f);
  if (col != 0) { fprintf(stderr, "ERROR: last line not terminated\n"); return 2; }
  if (row != n) { fprintf(stderr, "ERROR: %d rows of %d entries\n", row, n); return 2; }
  if (n > MAXN) { fprintf(stderr, "ERROR: too big\n"); return 2; }
  nw = (n + 63) / 64;
  long edges = 0;
  for (int i = 0; i < n; i++) {
    if (a[i][i]) { fprintf(stderr, "ERROR: nonzero diagonal\n"); return 1; }
    for (int j = 0; j < n; j++) {
      if (a[i][j] != a[j][i]) { fprintf(stderr, "ERROR: not symmetric\n"); return 1; }
      if (i < j && a[i][j]) edges++;
      int e = (i != j) && (ind ? !a[i][j] : a[i][j]);
      if (e) g[i][j >> 6] |= 1ULL << (j & 63);
    }
  }
  printf("file=%s\nn=%d edges=%ld symmetric=yes zero_diag=yes mode=%s T=%d threads=%d\n", argv[1], n, edges,
         ind ? "independent sets" : "cliques", T, omp_get_max_threads());
  fflush(stdout);
  unsigned long long tot[40]; memset(tot, 0, sizeof tot);
  double t0 = omp_get_wtime();
  #pragma omp parallel
  {
    unsigned long long cnt[40]; memset(cnt, 0, sizeof cnt);
    #pragma omp for schedule(dynamic, 1) nowait
    for (int v = 0; v < n; v++) {
      uint64_t C[MW]; int w = v >> 6; unsigned long long p = 0;
      for (int i = 0; i < nw; i++) C[i] = g[v][i];
      /* keep only neighbours with index > v */
      for (int i = 0; i < w; i++) C[i] = 0;
      C[w] &= (v & 63) == 63 ? 0ULL : ~((2ULL << (v & 63)) - 1);
      for (int i = w; i < nw; i++) p += __builtin_popcountll(C[i]);
      cnt[2] += p;
      if (p && 2 < T) rec(C, w, 1, cnt);
    }
    #pragma omp critical
    for (int i = 0; i < 40; i++) tot[i] += cnt[i];
  }
  tot[1] = n;
  for (int s = 1; s <= T; s++) printf("number of %s of size %d: %llu\n", ind ? "independent sets" : "cliques", s, tot[s]);
  printf("wall=%.1fs\n", omp_get_wtime() - t0);
  printf("RESULT: %s of size %d exists: %s\n", ind ? "independent set" : "clique", T, tot[T] ? "YES" : "no");
  return 0;
}
