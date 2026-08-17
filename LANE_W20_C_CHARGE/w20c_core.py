# w20c_core.py -- LANE W20_C (CHARGED SECTOR).  MACHINERY ONLY.  Prints nothing.
#
# Conventions inherited BYTE-IDENTICAL from LANE_W20_PRE / W-19:
#   Z_2 on each of L=12 links of tri_chain12.
#   G_v = prod_{l incident v} X_l ;  physical sector G_v = eta_v = (-1)^{q_v}.
#   H = -(1/g2) sum_{p in PLAQ} W_p  -  g2 sum_l X_l .
#   PLAQ = [1,2,3],[3,4,5],[7,8,9],[9,10,11],[0,1,4,6,7,10]   (pin_carrier BLOCK 7).
#
# REPRESENTATION.  Work in the X-eigenbasis:  X_l|u> = (-1)^{u_l}|u> ,  Z_l|u> = |u XOR e_l>.
# Then G_v|u> = (-1)^{|star(v) AND u|}|u>, so the physical sector for charge set q is
#      { u : d(u) = q }  =  u0 XOR Z_1 ,   dim 2^(L-V+1) = 32 .
# A gauge-invariant Pauli is  O_(a,c) = i^{|a AND c|} X^a Z^c  with c in Z_1 (Hermitian, square I).
#
# ALGEBRAIC ENTROPY.  For a Pauli algebra A with GF(2) symplectic space Vbar (after quotienting the
# operators that are SCALARS on the physical sector), Vbar = radical(dim r) + k hyperbolic pairs, so
# A = (+)_{2^r blocks} M_{2^k} and
#      S(rho|A) = H({p_s}) + sum_s p_s S(rho_s)
# with p_s = <psi|P_s|psi> and rho_s the 2^k x 2^k block matrix reconstructed from
#      t_{s,M} = <psi| P_s O_M |psi> .
# This is the standard algebraic (Casini-Huerta / Ohya-Petz) entropy; it reduces to the classical
# Shannon entropy when k = 0 and to the usual von Neumann entropy when r = 0.
import itertools, functools, math
import numpy as np

# ---------------------------------------------------------------- the carrier
V = 8
E = [(0, 7), (0, 1), (0, 2), (1, 2), (1, 3), (2, 3),
     (3, 4), (4, 5), (4, 6), (5, 6), (5, 7), (6, 7)]
L = len(E)
FULLMASK = (1 << L) - 1
star = [0] * V
for _i, (_a, _b) in enumerate(E):
    star[_a] |= 1 << _i
    star[_b] |= 1 << _i

def bits(m, n=L): return [i for i in range(n) if m >> i & 1]
def pop(m): return bin(m).count("1")

def dpart(m):
    r = 0
    for i in bits(m):
        a, b = E[i]
        r ^= (1 << a) ^ (1 << b)
    return r

CYC = sorted(m for m in range(1 << L) if dpart(m) == 0)          # 32 elements, dim 5
CIDX = {c: i for i, c in enumerate(CYC)}
NPHYS = len(CYC)
def delta(am): return functools.reduce(lambda x, v: x ^ star[v], bits(am, V), 0)
COC = sorted({delta(am) for am in range(1 << V)})                # 128 elements, dim 7
COCSET = set(COC)

S_LINKS = [1, 2, 3]; S_MASK = sum(1 << i for i in S_LINKS)
SIGMA   = [0, 4, 5]; SIG_MASK = sum(1 << i for i in SIGMA)
ENV     = [0, 4, 5, 6, 7, 8, 9, 10, 11]; ENV_MASK = sum(1 << i for i in ENV)
PLAQ    = [sum(1 << i for i in p) for p in
           ([1, 2, 3], [3, 4, 5], [7, 8, 9], [9, 10, 11], [0, 1, 4, 6, 7, 10])]
GRID    = [0.05, 0.10, 0.20, 0.30, 0.45, 0.60, 0.80, 1.00, 1.30, 1.70, 2.20, 3.00, 5.00]
DELTA_TOL = 0.10          # plateau tolerance, inherited
GATE      = 0.10          # emptiness gate, bits, pinned in pin_arms AXIS 6

# ---------------------------------------------------------------- GF(2) helpers
def gf2_basis(vs):
    b = []
    for v in vs:
        for x in b:
            v = min(v, v ^ x)
        if v:
            b = sorted(b + [v], reverse=True)
    return b

def gf2_rank(vs): return len(gf2_basis(vs))

def in_span(v, basis):
    for x in basis:
        v = min(v, v ^ x)
    return v == 0

def span_all(basis):
    out = [0]
    for b in basis:
        out += [x ^ b for x in out]
    return out

# ---------------------------------------------------------------- symplectic vectors
def sv(a, c): return (a << L) | c
def apart(v): return v >> L
def cpart(v): return v & FULLMASK
def sform(v, w):
    return (pop(apart(v) & cpart(w)) + pop(apart(w) & cpart(v))) & 1
def is_scalar(v):
    return cpart(v) == 0 and apart(v) in COCSET

SC_BASIS = gf2_basis([sv(star[v], 0) for v in range(V)])

# ---------------------------------------------------------------- Pauli operators with phase
def op_of(v):
    """canonical Hermitian representative: (a, c, phase) meaning phase * X^a Z^c, phase in {1,i,-1,-i}"""
    a, c = apart(v), cpart(v)
    return (a, c, (1j) ** (pop(a & c) % 4))

def op_mul(o1, o2):
    a1, c1, p1 = o1
    a2, c2, p2 = o2
    sgn = -1.0 if (pop(c1 & a2) & 1) else 1.0
    return (a1 ^ a2, c1 ^ c2, p1 * p2 * sgn)

OP_I = (0, 0, 1 + 0j)

# ---------------------------------------------------------------- physical sector
def u0_for(qmask):
    """minimum-weight u with d(u) = qmask; qmask must have even popcount."""
    best = None
    for u in range(1 << L):
        if dpart(u) == qmask:
            if best is None or pop(u) < pop(best):
                best = u
    return best

class Sector:
    """the 32-dim physical sector for a charge pattern; carries the u-table and the phase tables."""
    def __init__(self, charges):
        self.charges = tuple(sorted(charges))
        self.q = sum(1 << v for v in charges)
        assert pop(self.q) % 2 == 0, "odd charge is inadmissible: prod_v eta_v = +1"
        self.u0 = u0_for(self.q)
        self.U = np.array([self.u0 ^ z for z in CYC], dtype=np.int64)   # u of each basis index
        self.eta = np.array([1 - 2 * ((self.q >> v) & 1) for v in range(V)])
        self._pc = {}

    def parity(self, a):
        """(-1)^{|a AND u_z|} for every basis index z, as a +-1 float vector"""
        key = a
        if key not in self._pc:
            vals = np.array([1.0 - 2.0 * (pop(a & int(u)) & 1) for u in self.U])
            self._pc[key] = vals
        return self._pc[key]

    def perm(self, c):
        """index permutation induced by Z^c (c in Z_1): index of z XOR c"""
        return np.array([CIDX[z ^ c] for z in CYC], dtype=np.int64)

    def H(self, g2):
        M = np.zeros((NPHYS, NPHYS))
        for i, u in enumerate(self.U):
            M[i, i] = -g2 * (L - 2 * pop(int(u)))
        for p in PLAQ:
            pr = self.perm(p)
            for i in range(NPHYS):
                M[pr[i], i] += -1.0 / g2
        return M

    def ground(self, g2):
        w, vv = np.linalg.eigh(self.H(g2))
        return w[0], vv[:, 0].astype(complex), w

    def expect(self, op, psi):
        a, c, ph = op
        pr = self.perm(c)
        phase = self.parity(a)
        val = np.vdot(psi, ph * phase * psi[pr])
        return val

# ---------------------------------------------------------------- algebra structure (state-free)
class Alg:
    """Pauli algebra on the physical sector, from generator symplectic vectors."""
    _cache = {}

    def __init__(self, gens, name=""):
        self.name = name
        self.gens = tuple(sorted(set(gens)))
        B = gf2_basis(list(self.gens))
        self.rank = len(B)
        allv = span_all(B)
        scal = [v for v in allv if is_scalar(v)]
        K = gf2_basis(scal)
        self.dim_scalar = len(K)
        # representatives of Vbar: greedily extend K inside V
        rep = []
        cur = list(K)
        for v in allv:
            if not in_span(v, gf2_basis(cur)):
                rep.append(v)
                cur = gf2_basis(cur + [v])
        self.nbar = len(rep)
        assert self.nbar == self.rank - self.dim_scalar, (self.nbar, self.rank, self.dim_scalar)
        pairs, rad = symplectic_decompose(rep)
        self.pairs = pairs
        self.rad = rad
        self.k = len(pairs)
        self.r = len(rad)
        assert 2 * self.k + self.r == self.nbar, (self.k, self.r, self.nbar)
        self.dim = 1 << self.nbar            # linear dimension of the represented algebra
        self.maxent = self.k + self.r        # bits
        # central products C_T
        self.CT = []
        for T in range(1 << self.r):
            o = OP_I
            for i in bits(T, self.r):
                o = op_mul(o, op_of(self.rad[i]))
            self.CT.append((T, o))
        # block monomials
        self.MON = []
        for al in range(1 << self.k):
            for be in range(1 << self.k):
                o = OP_I
                for j in range(self.k):
                    if al >> j & 1:
                        o = op_mul(o, op_of(self.pairs[j][0]))
                    if be >> j & 1:
                        o = op_mul(o, op_of(self.pairs[j][1]))
                self.MON.append(((al, be), o))
        self._mats = None

    def block_mats(self):
        """abstract 2^k x 2^k matrices for each monomial, in the SAME generator order"""
        if self._mats is None:
            Z = np.array([[1, 0], [0, -1]], dtype=complex)
            X = np.array([[0, 1], [1, 0]], dtype=complex)
            I2 = np.eye(2, dtype=complex)
            out = {}
            for (al, be), _ in self.MON:
                M = np.array([[1.0 + 0j]])
                for j in range(self.k):
                    f = I2
                    if al >> j & 1: f = f @ Z
                    if be >> j & 1: f = f @ X
                    M = np.kron(M, f)
                out[(al, be)] = M
            self._mats = out
        return self._mats

    def state(self, sec, psi):
        """returns (p[s], rho_s) for each of the 2^r central sectors"""
        ev = {}
        for T, oT in self.CT:
            for M, oM in self.MON:
                o = op_mul(oT, oM)
                val = sec.expect(o, psi)
                assert abs(val.imag) < 1e-9, ("non-Hermitian expectation", val)
                ev[(T, M)] = val.real
        mats = self.block_mats()
        blocks = {}
        for s in range(1 << self.r):
            sgn = {}
            for T, _ in self.CT:
                sgn[T] = -1.0 if (pop(T & s) & 1) else 1.0
            R = np.zeros((1 << self.k, 1 << self.k), dtype=complex)
            for M, _ in self.MON:
                t = sum(sgn[T] * ev[(T, M)] for T, _ in self.CT) / (1 << self.r)
                R = R + t * mats[M]
            R = R / (1 << self.k)
            blocks[s] = R
        return blocks

    def entropy(self, sec, psi):
        blocks = self.state(sec, psi)
        tot = 0.0
        for s, R in blocks.items():
            w = np.linalg.eigvalsh((R + R.conj().T) / 2).real
            for x in w:
                if x > 1e-13:
                    tot += -x * math.log2(x)
                elif x < -1e-8:
                    raise AssertionError("negative eigenvalue %r in %s" % (x, self.name))
        return tot

def symplectic_decompose(basis):
    rest = [v for v in basis if v != 0]
    pairs = []
    while True:
        found = None
        for i in range(len(rest)):
            for j in range(i + 1, len(rest)):
                if sform(rest[i], rest[j]):
                    found = (i, j); break
            if found: break
        if not found: break
        i, j = found
        a, b = rest[i], rest[j]
        pairs.append((a, b))
        new = []
        for x_i, x in enumerate(rest):
            if x_i in (i, j): continue
            y = x
            if sform(x, b): y ^= a
            if sform(x, a): y ^= b
            new.append(y)
        rest = [y for y in new if y != 0]
    return pairs, gf2_basis(rest)

# ---------------------------------------------------------------- named algebras
def frag_gens(links):
    lm = sum(1 << l for l in links)
    g = [sv(1 << l, 0) for l in links]
    zs = gf2_basis([z for z in CYC if z and not (z & ~lm)])
    g += [sv(0, z) for z in zs]
    return g

A_FULL   = Alg([sv(1 << 1, 0), sv(1 << 2, 0), sv(1 << 3, 0), sv(0, S_MASK)], "FULL   A_S")
A_CENTRE = Alg([sv((1 << 1) | (1 << 2), 0), sv((1 << 1) | (1 << 3), 0)],      "CENTRE Z(A_S)")
A_BLOCK  = Alg([sv(1 << 1, 0), sv(0, S_MASK)],                                "BLOCK  {X_1,W_S}")
A_MAG    = Alg([sv(0, S_MASK)],                                               "MAG    {W_S}")
A_WIDE   = Alg(frag_gens([0, 1, 2, 3, 4, 5]),                                 "WIDE   S u Sigma")
CHANNELS = [("FULL", A_FULL), ("CENTRE", A_CENTRE), ("BLOCK", A_BLOCK), ("MAG", A_MAG)]

FRAG_P = [("F1", [0, 4, 5]), ("F2", [7, 8, 9]), ("F3", [11]), ("F4", [6, 10])]
FRAG_S = [("G1", [0, 4, 5]), ("G2", [6, 9]), ("G3", [7, 11]), ("G4", [8, 10])]
A_ENV  = Alg(frag_gens(ENV), "E_env")
A_SIG  = Alg(frag_gens(SIGMA), "Sigma")
AF_P   = [(n, Alg(frag_gens(f), n)) for n, f in FRAG_P]
AF_S   = [(n, Alg(frag_gens(f), n)) for n, f in FRAG_S]

_JOIN = {}
def join(a, b):
    key = (a.gens, b.gens)
    if key not in _JOIN:
        _JOIN[key] = Alg(list(a.gens) + list(b.gens), a.name + " v " + b.name)
    return _JOIN[key]

def MI(a, b, sec, psi):
    return a.entropy(sec, psi) + b.entropy(sec, psi) - join(a, b).entropy(sec, psi)

# ---------------------------------------------------------------- vacuity number D(F)
def Dnum(links):
    fm = sum(1 << l for l in links)
    return gf2_rank([c & S_MASK for c in COC if not (c & ~(S_MASK | fm))])

def cycdim(links):
    fm = sum(1 << l for l in links)
    return gf2_rank([z for z in CYC if not (z & ~fm)])

# ---------------------------------------------------------------- misc
def haar(sec, seed):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=NPHYS) + 1j * rng.normal(size=NPHYS)
    return v / np.linalg.norm(v)

def flux_sigma(sec):
    return int(np.prod([sec.eta[v] for v in (0, 1, 2)]))

def var_plaq(sec, psi):
    """Var_psi( sum_p W_p )"""
    tot = np.zeros_like(psi)
    for p in PLAQ:
        tot = tot + psi[sec.perm(p)]
    m1 = np.vdot(psi, tot).real
    m2 = np.vdot(tot, tot).real
    return m2 - m1 * m1, m1
