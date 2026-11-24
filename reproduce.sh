#!/bin/bash
# Rebuild the graph from Google's file, verify the hypotheses of Theorem 3, and check the graph with exact programs.
# Needs: python3 (standard library), gcc (with OpenMP for check_enumerate and check_branch_bound_fable).
# Default run: about 15 to 30 minutes depending on machine load (the Russian-doll clique search is the longest part, 12 to 24 minutes, one core).
# What this script can and cannot show: it proves that the graph has no K9 (three exact programs), that the hypotheses of
# Theorem 3 hold, and that the G block has no K9 and no independent 9-set. It does NOT search the finished graph for an
# independent 17-set; that side rests on Theorem 3. The optional long mode only starts such a search (see below).
#   bash reproduce.sh long     (adds a direct search for an independent 17-set; no run time is known, it was never completed)
set -e
cd "$(dirname "$0")"
sha256sum -c SHA256SUMS
echo "== Rebuild"
python3 build.py theorem3_R9_17_1413.txt input_google_R9_11_617.txt rebuilt_1413.txt
cmp rebuilt_1413.txt R9_17_1413.txt
echo "rebuilt graph is byte-identical to R9_17_1413.txt"
echo "== Hypotheses of Theorem 3 (about 1 minute)"
python3 verify_theorem3.py theorem3_R9_17_1413.txt input_google_R9_11_617.txt R9_17_1413.txt
gcc -O3 -o check_russian_doll check_russian_doll.c
gcc -O3 -fopenmp -o check_enumerate check_enumerate.c
gcc -O3 -fopenmp -o check_branch_bound_fable check_branch_bound_fable.c
echo "== G block (rows and columns 0..563, Mathon's graph): no K9 and no independent 9-set"
python3 - <<'PY'
import re
d = [int(c) for c in re.sub(r'[^01]', '', open("R9_17_1413.txt").read())]
N = int(round(len(d) ** 0.5))
open("G564.txt", "w").write("".join(" ".join(map(str, d[i * N:i * N + 564])) + "\n" for i in range(564)))
PY
./check_enumerate G564.txt 9 c | tee out_G564_c.txt; grep -q "exists: no" out_G564_c.txt
./check_enumerate G564.txt 9 i | tee out_G564_i.txt; grep -q "exists: no" out_G564_i.txt
echo "== Clique side, method 1: Russian-doll search (the line 'clique number = N' is the exact clique number)"
./check_russian_doll R9_17_1413.txt 9 c | tee out_1413_c.txt; grep -q "exists: no" out_1413_c.txt
echo "== Clique side, method 2: plain enumeration with no pruning, counts every size"
./check_enumerate R9_17_1413.txt 9 c | tee out_1413_ce.txt; grep -q "exists: no" out_1413_ce.txt
echo "== Clique side, method 3: branch and bound with colour classes"
./check_branch_bound_fable R9_17_1413.txt 9 c 4 | tee out_1413_bb.txt; grep -q "exists: no" out_1413_bb.txt
echo "== Positive controls: a K8 exists, an independent 16-set exists"
./check_branch_bound_fable R9_17_1413.txt 8 c 4 | tee out_1413_bb8.txt; grep -q "exists: YES" out_1413_bb8.txt
./check_branch_bound_fable R9_17_1413.txt 16 i 4 | tee out_1413_bb16i.txt; grep -q "exists: YES" out_1413_bb16i.txt
if [ "$1" = "long" ]; then
  echo "== Independent-set side, direct search for an independent 17-set (long; never completed by us)"
  ./check_branch_bound_fable R9_17_1413.txt 17 i 8 | tee out_1413_bb17i.txt; grep -q "exists: no" out_1413_bb17i.txt
fi
echo "ALL CHECKS PASSED"
