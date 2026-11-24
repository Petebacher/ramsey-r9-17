"""Check the hypotheses of Theorem 3 (Xu, Shao, Radziszowski, SIAM J. Discrete Math. 25 (2011) 394-400) for a graph F,
standard library only. Usage: python3 verify_theorem3.py theorem3_R9_17_1413.txt input_google_R9_11_617.txt R9_17_1413.txt
Checks, directly on the finished graph F:
  1. F is a symmetric 0/1 matrix with zero diagonal, n(F) = n(G) + n(H) + m.
  2. The G block equals the stated rows of the Google file (and equals Paley(n(G)) when g_source is paley, and equals the
     Taylor graph of Paley(q) under the stated vertex map when a taylor_map is given); the H block equals the stated rows.
  3. The only G-H edges are the matching u_i -- v_i, and each new vertex z_i is joined to exactly N_G(u_i) and N_H(v_i).
  4. M = G[U] equals H[V] under u_i -> v_i (a common induced subgraph on m vertices).
  5. The new vertices induce M[W1] + M[W2] with no W1-W2 edge, both W1 and W2 nonempty.
  6. Split condition: for i = 1, 2 and every v of M outside Wi, M[Wi + v] has no K_(k-1).
  7. H is a (k,t)-graph: no K_k and no independent t-set (exact, only when H is small; for G = H = the 564-vertex block it is
     checked that H is the same graph as G, and the (k,t) property is left to the C programs).
  8. Counting: n(F) = n(G) + n(H) + n(M), and the theorem then gives R(k, s+t-1) >= n(F) + 1.
NOT checked here (too slow in pure Python; done by the C programs, see README): G has no K_k and no independent s-set.
Exit status 0 only if every check above passes."""
import hashlib, os, re, sys


def matrix(path):
    d = [int(c) for c in re.sub(r'[^01]', '', open(path).read())]
    n = int(round(len(d) ** 0.5))
    assert n * n == len(d), 'not square'
    return [d[i * n:(i + 1) * n] for i in range(n)]


def spec(path):
    sp = {}
    for line in open(path):
        line = line.split('#')[0].split()
        if line:
            sp[line[0]] = line[1:]
    return sp


def omega(adj, verts):
    """exact clique number of the graph on `verts`; adj[v] is a bitmask of the neighbours of v"""
    best = 0

    def rec(size, cand):
        nonlocal best
        if size > best:
            best = size
        while cand:
            if size + bin(cand).count('1') <= best:
                return
            v = cand.bit_length() - 1
            cand &= ~(1 << v)
            rec(size + 1, cand & adj[v])
    m = 0
    for v in verts:
        m |= 1 << v
    rec(0, m)
    return best


def mask(rows):
    return [sum(1 << j for j, x in enumerate(r) if x) for r in rows]


sp = spec(sys.argv[1])
A = matrix(sys.argv[2])
F = matrix(sys.argv[3])
k, s, t = int(sp['k'][0]), int(sp['s'][0]), int(sp['t'][0])
nG, nH, h0 = int(sp['nG'][0]), int(sp['nH'][0]), int(sp['h_start'][0])
pairs = [tuple(map(int, p.split(':'))) for p in sp['pairs']]
W1 = set(map(int, sp['W1']))
m = len(pairs)
W2 = set(range(m)) - W1
U = [u for u, v in pairs]
V = [v for u, v in pairs]
N = len(F)
ok = True


def check(name, cond):
    global ok
    print(('PASS  ' if cond else 'FAIL  ') + name)
    ok = ok and bool(cond)


check('1. F is %d x %d, symmetric, zero diagonal; n(F) = %d + %d + %d' % (N, N, nG, nH, m),
      N == nG + nH + m and all(F[i][j] == F[j][i] for i in range(N) for j in range(N))
      and all(F[i][i] == 0 for i in range(N)))
G = [r[:nG] for r in F[:nG]]
H = [r[nG:nG + nH] for r in F[nG:nG + nH]]
check('2a. G block = rows and columns 0..%d of the Google file' % (nG - 1), G == [r[:nG] for r in A[:nG]])
if sp['g_source'][0] == 'paley':
    sq = {(x * x) % nG for x in range(1, nG)}
    prime = all(nG % q for q in range(2, int(nG ** .5) + 1))
    check('2b. G block = Paley(%d): i ~ j iff i - j is a nonzero square mod %d (%d is prime and 1 mod 4)' % (nG, nG, nG),
          prime and nG % 4 == 1 and all(G[i][j] == (1 if (i - j) % nG in sq else 0) for i in range(nG) for j in range(nG)))
if 'taylor_map' in sp:
    q = int(sp['taylor_q'][0])
    f = {}
    for line in open(sp['taylor_map'][0]):
        if line.strip():
            i, j = map(int, line.split()); f[i] = j
    sqs = {(x * x) % q for x in range(1, q)}
    prime_q = all(q % r for r in range(2, int(q ** .5) + 1)) and q % 4 == 1

    def tay(i, j):   # Taylor(Paley(q)): 0..q-1 = x, q = apex, q+1..2q = x', 2q+1 = apex'
        if i == j:
            return 0
        if i > j:
            i, j = j, i
        if j == q:
            return 1 if i < q else 0
        if j == 2 * q + 1:
            return 1 if i > q else 0
        if i < q and j < q:
            return 1 if (j - i) % q in sqs else 0
        if i == q:   # apex versus a primed vertex
            return 0
        if i < q:   # x ~ y' iff y - x is a non-square
            d = (j - q - 1 - i) % q
            return 0 if (d == 0 or d in sqs) else 1
        return 1 if (j - i) % q in sqs else 0       # both primed
    nT = 2 * q + 2
    check('2c. G block = Taylor graph of Paley(%d) (built from scratch here) under the explicit vertex map in %s (%d is prime and 1 mod 4)' % (q, sp['taylor_map'][0], q),
          prime_q and nT == nG and sorted(f) == list(range(nT)) and sorted(f.values()) == list(range(nT))
          and all(G[f[i]][f[j]] == tay(i, j) for i in range(nT) for j in range(nT)))
check('2d. H block = rows and columns %d..%d of the Google file' % (h0, h0 + nH - 1),
      H == [r[h0:h0 + nH] for r in A[h0:h0 + nH]])
cross = {(i, j - nG) for i in range(nG) for j in range(nG, nG + nH) if F[i][j]}
check('3a. the G-H edges are exactly the %d matching edges u_i -- v_i' % m,
      cross == set(pairs) and len(set(U)) == m and len(set(V)) == m)
z0 = nG + nH
twins = all(F[z0 + i][:nG] == G[U[i]] and F[z0 + i][nG:nG + nH] == H[V[i]] for i in range(m))
check('3b. each z_i has exactly the neighbours N_G(u_i) in G and N_H(v_i) in H (so it is not joined to u_i or v_i)', twins)
M = [[G[U[i]][U[j]] for j in range(m)] for i in range(m)]
check('4. G[U] = H[V] entry by entry under u_i -> v_i; M has %d vertices and %d edges' % (m, sum(map(sum, M)) // 2),
      all(G[U[i]][U[j]] == H[V[i]][V[j]] for i in range(m) for j in range(m)))
Z = [[F[z0 + i][z0 + j] for j in range(m)] for i in range(m)]
check('5. Z = M[W1] + M[W2] with no W1-W2 edge; |W1| = %d, |W2| = %d, both nonempty' % (len(W1), len(W2)),
      bool(W1) and bool(W2)
      and all(Z[i][j] == (M[i][j] if (i in W1) == (j in W1) else 0) for i in range(m) for j in range(m)))
Mm = mask(M)
if os.path.exists('common_part_split_R9_17_1413.txt'):
    plain = [tuple(map(int, l.split())) for l in open('common_part_split_R9_17_1413.txt') if l.strip() and not l.startswith('#')]
    check('4b. common_part_split_R9_17_1413.txt (vertex, side) agrees with the pairs and split in the data file',
          plain == [(U[i], 1 if i in W1 else 2) for i in range(m)])
oW = {name: omega(Mm, sorted(Wi)) for name, Wi in (('W1', W1), ('W2', W2))}
check('5b. clique number of M[W1] and M[W2]: %d and %d (the paper\'s simpler sufficient form needs at most k-3 = %d)'
      % (oW['W1'], oW['W2'], k - 3), max(oW.values()) <= k - 3)
worst = {}
for name, Wi in (('W1', W1), ('W2', W2)):
    worst[name] = max(omega(Mm, sorted(Wi | {v})) for v in range(m) if v not in Wi)
check('6. for every v outside Wi, M[Wi + v] has no K_%d: largest clique found %d (W1) and %d (W2), need at most %d'
      % (k - 1, worst['W1'], worst['W2'], k - 2), max(worst.values()) <= k - 2)
print('      (clique number of M itself: %d; survey item 2.3.h would need at most %d)' % (omega(Mm, range(m)), k - 2))
if nH <= 100:
    Hm = mask(H)
    Hc = mask([[0 if i == j else 1 - H[i][j] for j in range(nH)] for i in range(nH)])
    wH = omega(Hm, range(nH))
    aH = omega(Hc, range(nH))
    check('7. H is a (%d,%d)-graph: clique number %d (need < %d), independence number %d (need < %d)' % (k, t, wH, k, aH, t),
          wH < k and aH < t)
else:
    check('7. H is the same graph as G (identical block, G = H, t = s = %d); its (%d,%d) property is the C-program check of G' % (t, k, t),
          G == H and s == t)
cf = ''.join(''.join(map(str, r)) + '\n' for r in F)
print('      sha256 of F with each row written without spaces: ' + hashlib.sha256(cf.encode()).hexdigest())
check('8. counting: n(G) + n(H) + n(M) = %d + %d + %d = %d = n(F); s + t - 1 = %d + %d - 1 = %d'
      % (nG, nH, m, nG + nH + m, s, t, s + t - 1), nG + nH + m == N)
print('NOT CHECKED HERE: G (= H) has no K_%d and no independent %d-set (see README, checks table)' % (k, s))
if ok:
    print('ALL HYPOTHESES VERIFIED. If G is a (%d,%d)-graph, Theorem 3 gives R(%d,%d) >= %d.' % (k, s, k, s + t - 1, N + 1))
sys.exit(0 if ok else 1)
