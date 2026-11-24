# R(9,17) ≥ 1413

This repository gives a graph on 1412 vertices that contains no clique of 9 vertices and, by a published theorem, no independent set of 17 vertices. It shows that the Ramsey number R(9,17) is at least 1413.

**Please read this first: what is checked by computer and what is not.** The absence of a clique of 9 vertices is checked directly on the finished graph by three independently written exact programs. The absence of an independent set of 17 vertices is **not** checked directly: no program has searched the finished graph for such a set, and a full search at this size is out of reach for us. That side rests on Theorem 3 of Xu, Shao and Radziszowski (2011), whose hypotheses we checked by machine on the finished graph. Details are in the Verification section.

Companion results: R(9,11) ≥ 620 at https://github.com/Petebacher/ramsey-r9-11 and R(10,12) ≥ 857 at https://github.com/Petebacher/ramsey-r10-12.

## Results

| Bound | Graph file | Vertices | Previous best public bound we could find | Gain |
|---|---|---|---|---|
| R(9,17) ≥ 1413 | `R9_17_1413.txt` | 1412 | 1411 (Xu, Xie, Radziszowski 2004, as listed in the survey) | +2 |

Prior art. The dynamic survey *Small Ramsey Numbers* (S. P. Radziszowski, DS1 revision 18, April 2026), note 2.2(g), states R(9,17) ≥ 1411 and cites Xu Xiaodong, Xie Zheng and S. P. Radziszowski, "A constructive approach for the lower bounds on the Ramsey numbers R(s,t)", J. Graph Theory 47 (2004) 231-239 (their Theorem 4 with p = 9, which glues two copies of a 564-vertex (9,9)-graph along the 282-vertex neighbourhood of one vertex; for the graph below, which is regular of degree 281, that means its complement). The value has been listed in the survey since at least the 2006 revision. We found nothing higher in the survey, the Leaps in Bounds tracker, Jang's gluing repository, the Quantavane repository, Google Research's `ramsey_number_bounds` folder, arXiv, or a general web search. These are the best public bounds we could find as of 24 November 2026. We searched zbMATH Open (24 November 2026: Ramsey lower bounds, Paley and Taylor graph constructions, R(9,17)) and found nothing above 1411; the closest later work, Exoo and Tatarevic, "New lower bounds for 28 classical Ramsey numbers" (2015), does not include R(9,17), and DS1 revision 18 still lists 1411. We did not search MathSciNet (subscription only), social media, MathOverflow or mailing lists, so this is not a claim about everything that is known.

## Where the graph comes from

The graph uses **Theorem 3** of Xu Xiaodong, Shao Zehui and S. P. Radziszowski, "More constructive lower bounds on classical Ramsey numbers", SIAM J. Discrete Math. 25 (2011) 394-400. We read the theorem in the authors' preprint (dated 14 September 2010); the theorem number is as in the preprint, and we could not check it against the journal version. The theorem says:

> Let k ≥ 4 and s, t ≥ 2. Let G be a (k,s)-graph (no K_k, no independent s-set) and H a (k,t)-graph. Let M be a graph isomorphic to an induced subgraph of both G and H. If the vertex set of M can be split into two nonempty sets W1 and W2 so that for i = 1, 2 and every vertex v of M outside Wi the graph M[Wi + v] has no K_(k−1), then R(k, s+t−1) ≥ n(G) + n(H) + n(M) + 1.

(The authors also note a simpler sufficient form: it is enough that neither M[W1] nor M[W2] contains a K_(k−2). Our split satisfies both forms.)

The inputs are not new. G and H are both Mathon's 564-vertex (9,9)-graph (A. Mathon, 1987, giving R(9,9) ≥ 565), which is rows 0 to 563 of Google's 616-vertex R(9,11) file (the identity labelling; G and H are two copies of the same graph, which Theorem 3 allows). This graph is the Taylor graph of the Paley graph on 281 vertices: two copies of Paley(281) (281 is prime and 1 mod 4; vertices i and j are adjacent when i − j is a nonzero square mod 281), each with an extra vertex joined to all of its copy, and each vertex x joined to the vertices y' of the other copy for which y − x is a non-square. It is 281-regular; two adjacent vertices have exactly 140 common neighbours; every vertex has exactly one "antipode", the unique non-neighbour that shares no neighbour with it. We checked that the file's block equals this graph exactly under an explicit vertex map (`taylor_map_G564.txt`, checked by `build.py` and `verify_theorem3.py`).

**The new part is the common part M (284 vertices) and its split.** Take an edge ab = {0, 1} of G and its antipodal edge a'b' = {281, 282} (the antipodes of 0 and 1; it is an edge). Let C(xy) be the set of common neighbours of x and y. Then

- M = (N[a] ∩ N[b]) ∪ (N[a'] ∩ N[b']), where N[x] is the closed neighbourhood; it has 142 + 142 = 284 vertices;
- W1 = {a, b} ∪ C(a'b') and W2 = C(ab) ∪ {a', b'}, 142 vertices each;
- there is no edge between {a, b} and C(a'b'), and none between {a', b'} and C(ab). The new block of the finished graph (M[W1] plus M[W2], with the M edges between the two sides removed) has components of sizes 2, 140, 140 and 2, which forces the split.

| | R(9,17) ≥ 1413 |
|---|---|
| k, s, t | 9, 9, 9 |
| G, H | Mathon's 564-vertex graph = Taylor(Paley(281)), a (9,9)-graph, used twice |
| M | 284 vertices, 20,022 edges, clique number 8 |
| split | 142 + 142 |
| clique number of M[W1] and of M[W2] | 6 (the simpler form needs at most k−3 = 6) |
| largest clique in any M[Wi + v] | 7 (limit k−2 = 7; the theorem forbids K8) |
| count | 564 + 564 + 284 = 1412 |
| s + t − 1 | 9 + 9 − 1 = 17 |

M itself contains a K8, so the simpler survey rule (item 2.3.h, a common part without K_(k−1)) does not apply here; the split in Theorem 3 is what makes it work, with nothing to spare (7 is exactly the limit). The earlier bound 1411 also glues two copies of a 564-vertex (9,9)-graph (the complement of this one) with a 282-vertex common part; the gain here is that the common part has 284 vertices.

The pairing and the split are plain files: `theorem3_R9_17_1413.txt` (the format read by the scripts: 284 pairs u:v, here always u:u because G = H, and the split) and `common_part_split_R9_17_1413.txt` (one line per vertex of the common part: the vertex of G and its side, 1 or 2). The verifier checks that the two files agree.

## The graph F built from the data

F has the vertices of G (0 to 563), then the vertices of H (564 to 1127), then 284 new vertices z_1 ... z_284 (1128 to 1411), one per pair (u_i, v_i):
- the edges of G and of H are kept;
- u_i is joined to v_i (a matching; these are the only edges between G and H);
- z_i is joined to the neighbours of u_i in G and to the neighbours of v_i in H (not to u_i or v_i themselves);
- z_i and z_j are joined when u_i u_j is an edge of G and i, j are on the same side of the split.
`build.py` does this using the Python standard library only. It also builds the Taylor graph of Paley(281) from scratch and compares it with the file's block under the vertex map.

## Why it works (Theorem 3, short form)

- **No K_k.** A clique that meets both G and H has exactly one vertex u_j in G and one v_j in H, since only matching edges join them. Its other k−2 vertices cannot be in both parts of the split, so they are new vertices z_i with all their indices i in one Wi, and those indices together with j form a K_(k−1) in M[Wi + j] (or in M[Wi], when j itself is in Wi). The split condition excludes this. A clique missing H or G is mapped to a clique of G (or H) of the same size by replacing each z_j by u_j (v_j), so it has fewer than k vertices.
- **No independent set of size s + t − 1.** Split an independent set I into its parts in G, in H, in W1 and in W2, and sort the G and H vertices by whether their partner z is in I. The parts that sit next to G embed in G as one independent set, giving at most s − 1 vertices; the others embed in H, giving at most t − 1. So |I| ≤ s + t − 2.

The full proof is in the paper. The clique side is also checked directly (below); **the independent-set side is not**, and for it you rely on the proof. The proof was independently re-derived by a separate AI review, which found that it covers this case (G = H, the identity map, and 9 + 9 − 1 = 17 are all allowed).

## Verification

Everything below was run on the files in this repository (the sha256 values are in `SHA256SUMS`).

**File format.** `R9_17_1413.txt` is 1412 lines, each a row of 1412 entries 0/1 separated by single spaces. There is no header line. All checkers read every `0` and `1` character of the file and nothing else, so do not add a header line. The same matrix with each row written without spaces (one string of 1412 digits per line) has sha256 `a3ceb36128af8fd72bda92978292b24736299d4fce48b1c175c67beeec9c6901`; `verify_theorem3.py` prints this value.

**1. Hypotheses of Theorem 3** (`verify_theorem3.py`, standard library only, about one minute): F is symmetric with zero diagonal; its G block is exactly the stated rows of Google's file and equals Taylor(Paley(281)) under the vertex map; its H block is the same rows; the G-H edges are exactly the matching; every z_i has exactly the neighbourhoods N_G(u_i), N_H(v_i); M is the same graph in G and in H; Z = M[W1] + M[W2] with no W1-W2 edge; the clique numbers of M[W1] and M[W2] are 6; for every v outside Wi the largest clique of M[Wi + v] is 7 (an exact search, for all 284 choices of v), so none has 8 vertices; the counts add up. All pass (`logs/verify_theorem3_R9_17_1413.log`). One further, independent run of the same hypotheses by a separate AI review, with its own code, gave the same numbers (`logs/fable_verify_c1.log`; that script is not shipped).

**2. G is a (9,9)-graph.** Done with `check_enumerate` (exhaustive enumeration, no pruning): the 564-vertex block has no K9 (210,783,720 cliques of size 8, none of size 9) and no independent 9-set (same counts). Logs `logs/enumerate_G564_c9.log`, `logs/enumerate_G564_i9.log` (from the R(9,11) release) and the rerun `logs/check_enumerate_G564_*.log`.

**3. Direct exact checks of the finished graph for cliques.** Every run says "no K9":

| Program | Method | R9_17_1413.txt (1412 vertices) |
|---|---|---|
| `check_russian_doll.c` | Östergård's Russian-doll method, exact | clique number exactly 8 (15,472,846,906 nodes; 728 s in the first run, 1426 s cpu in the rerun on a busy machine; same node count) |
| `check_enumerate.c` | plain enumeration with no pruning at all, counts every size | 4,461,521,120 cliques of size 8 and 0 of size 9 (207 s wall, 8 threads) |
| `check_branch_bound_fable.c` | branch and bound with colour classes (a third, separately written method) | no K9 (812,097,401 nodes in the first run; the same node count in the rerun, 111 s wall on 4 threads) |

Counts from `check_enumerate.c`, which anyone can match against a rerun: 1,412 vertices; 328,038 edges; 21,420,284 triangles; 493,640,140 cliques of size 4; 4,267,061,120 of size 5; 13,668,755,520 of size 6; 14,812,369,600 of size 7; 4,461,521,120 of size 8; 0 of size 9.

**4. Positive controls.** The searches find sets one size below the target: the graph has cliques of size 8, and `check_branch_bound_fable` finds and verifies an explicit independent 16-set (so the independence number is at least 16; if the theorem holds it is exactly 16). Logs `logs/check_branch_bound_fable_R9_17_1413_control_*.log`.

**5. Independent sets of size 17: not checked directly.** There is no program run on the finished graph that shows it has no independent set of 17 vertices. Searching 1412 vertices for such a set exactly is far beyond what we ran (an exact independence-number search on the 856-vertex graph of our R(10,12) release took 4.2 hours; this graph has 1412 vertices and independence number at least 16). The claim for this side is the proof of Theorem 3 together with the machine-checked hypotheses in item 1 and the check in item 2 that G is a (9,9)-graph. The small-scale instances of the same recipe on Taylor(Paley(13)) and Taylor(Paley(17)) (72 and 92 vertices) were checked exactly on both sides by the separate review and sit exactly at the theorem's limits (clique number 4, independence number 8); that code is not shipped. `reproduce.sh long` starts a direct search for an independent 17-set; we have never run it to completion and give no run time.

How independent the checks are: the three clique programs use three different methods (Russian doll; unpruned enumeration; colour-class branch and bound), and they come from three separate pieces of work: the Russian-doll program was written by the session that built our graphs, the enumeration program by the AI review of our R(9,11) release, and the branch-and-bound program by the AI review of this graph.

**Logs.** `logs/` holds the finished runs. Notes:
- The logs named `rdoll_*`, `enumerate_F1412_c9.log`, `branch_bound_fable_*` and `fable_verify_c1.log` were made on a file with identical content (the same matrix, rows written without spaces, sha256 above) under the working name `F1412.txt`; the first line of some logs shows a working name. In branch_bound_fable_F1412_c9.log a directory path was removed from that line; nothing else in it was changed. Those runs used builds of the same programs that differ from the shipped ones only in the array-size constant MAXN and, for the first Russian-doll log, in the printing labels: `rdoll_main_F1412_clique9.log` used the Russian-doll program before the printing fix described next, with MAXN 1600; `rdoll_decompose_F1412_c9.log` and `enumerate_F1412_c9.log` used the corrected programs with MAXN 1472. The shipped programs have MAXN 1600. (The sister repositories' copies of these two programs have MAXN 1100 and 1024, which are too small for 1412 vertices; we raised it here and say so.)
- `rdoll_main_F1412_clique9.log` was written by the version of the Russian-doll program before a printing fix. In clique mode that version called the result "independence number" (it works on the complement, so the number it prints is the clique number), and it printed "K4 count: 0" even when the count was skipped for n > 300. The shipped `check_russian_doll.c` fixes only these printed labels; the search code is unchanged. The change is in `logs/check_russian_doll_label_fix.diff` (made for the R(10,12) release). Read "independence number" in that log as "clique number", and ignore "K4 count: 0".
- The logs whose names start with `check_` were produced by the shipped programs on the shipped file on 24 November 2026. They show the corrected labels.

`reproduce.sh` runs the checksum check, the rebuild with byte comparison, the hypothesis verifier, the G-block checks, the three clique-side checks and the positive controls (about 15 to 30 minutes depending on machine load; 27 minutes on a busy machine). `logs/reproduce_default_run.log` is the output of one full default run on 24 November 2026 (all checks passed). `reproduce.sh long` additionally starts the direct independent-17 search described in item 5.

## Files

- `R9_17_1413.txt`: the graph (format above).
- `input_google_R9_11_617.txt`: Google's 616-vertex R(9,11) graph, redistributed unchanged (renamed from `R(9, 11) >= 617.txt`) with attribution; see License. Rows 0 to 563 are the block used.
- `theorem3_R9_17_1413.txt`: pairing and split. `common_part_split_R9_17_1413.txt`: the same split as a plain list (vertex, side). `taylor_map_G564.txt`: the explicit vertex map from Taylor(Paley(281)) to the Google block.
- `build.py`: rebuilds the graph and checks the Taylor description. `verify_theorem3.py`: the hypothesis checker.
- `check_russian_doll.c` (usage `file T i|c`), `check_enumerate.c` (usage `file T c|i`, needs OpenMP), `check_branch_bound_fable.c` (usage `file T c|i [threads]`, needs OpenMP). The first two have MAXN = 1600.
- `reproduce.sh`, `SHA256SUMS`, `logs/`.

## Credit and method

The starting graph is Mathon's 564-vertex graph, as published in Google Research's R(9,11) file. Theorem 3 is due to Xu, Shao and Radziszowski, and the earlier bound 1411 to Xu, Xie and Radziszowski. The new part is the 284-vertex common part and its split, and the checks.

**How we found the common part.** We found it by taking apart Google Research's R(7,13) ≥ 513 graph (google-research/google-research, `ramsey_number_bounds`, 22 October 2026). That graph consists of two copies of a 204-vertex graph H204, which is the Taylor graph of Paley(101), and 104 further vertices; its common part is exactly the same two-antipodal-edge form (104 = 2 × (λ + 2) vertices, with λ = 50 common neighbours per edge), with the same two sides. Rebuilding Google's graph from this recipe reproduced the file with no differing entries, in our own check and again in a separate AI review's check with its own code; neither rebuild is part of this repository, and anyone can repeat it from Google's file. Applying the same recipe to the Taylor graph of Paley(281) gives 2 × (140 + 2) = 284. We did not discover the recipe independently of Google's graph. The general statement is that for every prime power q = 1 mod 4 with q ≥ 9 the two-antipodal-edge common part of Taylor(Paley(q)) satisfies the hypotheses of Theorem 3 automatically and gives R(w+2, 2w+3) ≥ 5q + 8, where w is the clique number of Paley(q); we checked it exactly only on small cases, and here only for q = 281. This is exactly 2 more than the value 5q + 6 that Theorem 4 of Xu, Xie and Radziszowski (2004) gives for these graphs (the Mathon-Shearer graphs). It does not improve Theorem 4 as a general statement, which holds for every (p,p)-graph.

This work was directed by Peter J. Chrisbacher. The search, construction and checks were carried out with AI coding agents (Anthropic's Claude Code), working under fixed rules: every claim needs independently written exact checks and a check of published results before any announcement. A separate AI review re-derived the theorem's proof, re-checked every hypothesis with its own code, found the prior art above and wrote the branch-and-bound program shipped here as `check_branch_bound_fable.c`. The enumeration and Russian-doll programs are shared with our R(9,11) and R(10,12) releases.

## Changelog

- **v1.0 (24 November 2026).** First release: R(9,17) ≥ 1413.

## Glossary

- **Ramsey number R(s,t):** the smallest n such that every graph on n vertices has a clique of size s or an independent set of size t. A graph on n vertices with no K_s and no independent t-set shows R(s,t) ≥ n + 1.
- **Clique, K_k:** a set of k vertices that are all joined to each other. **Independent set:** a set of vertices with no edges among them.
- **(k,s)-graph:** a graph with no K_k and no independent s-set.
- **Closed neighbourhood N[x]:** x together with its neighbours. **C(xy):** the common neighbours of x and y.
- **Induced subgraph M[W]:** the graph on the vertices W with all edges of M among them.
- **Taylor graph of Paley(q):** see above; it has 2(q+1) vertices. **Antipode:** in a Taylor graph, the unique non-neighbour of a vertex that has no common neighbour with it.
- **Twin:** a vertex with the same neighbours as another vertex (the new vertices z_i are such copies of u_i and v_i, with some edges removed by the split).

## License

- Code (`build.py`, `verify_theorem3.py`, the checkers, `reproduce.sh`): MIT.
- The new graph, the pairing, split and vertex-map files and this text: CC BY 4.0.
- `input_google_R9_11_617.txt` is Google's file, redistributed unchanged with attribution under the terms of google-research/google-research (Apache 2.0 for source files, CC BY 4.0 for datasets; the Apache text is in `LICENSES/`). The graph `R9_17_1413.txt` contains Google's 564-vertex block twice, unchanged (plus new vertices and edges).

## How to cite

See `CITATION.cff`, or cite the Zenodo record for this repository once it is published.
