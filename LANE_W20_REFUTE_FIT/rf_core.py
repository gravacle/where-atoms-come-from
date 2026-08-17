# rf_core.py -- LANE W20_REFUTE_FIT.  INDEPENDENT MACHINERY.
#
# WHY THIS FILE EXISTS AT ALL.  I am the refuter of the FITTING charge for W-20 attempt two.
# I do not import LANE_W20_R_LEDGER/core_w20r.py or LANE_W20_C_CHARGE/w20c_core.py.  If I
# reproduce their numbers with a differently-built instrument, the numbers are not an artefact
# of one lane's algebra code; if I fail to, that is itself the finding.  The ONLY things taken
# from the other lanes are the CARRIER, the PLAQUETTE SET, the HAMILTONIAN CONVENTION, the
# REGION/SURFACE and the COUPLING GRID -- because a refuter that changes the carrier is not
# refuting anything, it is running a different experiment.  Those are listed as INHERITED below.
#
# INHERITED VERBATIM (declared, not re-derived):
#   carrier tri_chain12: V=8, L=12, cubic, girth 3, EDGES below (W20_PRE pin_carrier).
#   PLAQ = [1,2,3],[3,4,5],[7,8,9],[9,10,11],[0,1,4,6,7,10]   (W20_PRE BLOCK 7)
#   H = -(1/g2) sum_p W_p - g2 sum_l X_l                       (W-19 convention)
#   region A = {0,1,2}, S = links{1,2,3}, Sigma = delta(A) = links{0,4,5}
#   13-point grid [0.05 ... 5.00]                              (W20_PRE BLOCK 5)
#   Z_2 conventions: Z|j> = (-1)^j|j>, X|j> = |j+1>, G_v = prod_{l ~ v} X_l.
#
# BUILT INDEPENDENTLY HERE (different method from both lanes):
#   * the physical sector is built as an EXPLICIT 4096 x 32 ISOMETRY U, not as a coset backend.
#     Every operator is then a dense 32x32 matrix U^dag O U.  No Pauli-label bookkeeping, no
#     effective-label quotient, no symplectic decomposition.  Dense linear algebra only.
#   * the algebra entropy of A_S is NOT computed from a GF(2) symplectic normal form.  It is
#     computed from the CONCRETE central projectors and a CONCRETE Bloch vector:
#         Z(A_S) = alg{ c1 = X_1 X_2 , c2 = X_1 X_3 }        (abelian, 4 minimal projectors)
#         block  = alg{ bx = X_1 , bz = W_S = Z_1 Z_2 Z_3 }  (bx, bz anticommute -> one qubit)
#         P_k = (1 + s1 c1)/2 (1 + s2 c2)/2,  p_k = <psi|P_k|psi>
#         rho_k = (I + x_k sx + y_k sy + z_k sz)/2 with
#           x_k = <P_k bx>/p_k , z_k = <P_k bz>/p_k , y_k = <P_k (i bx bz)>/p_k
#         S(rho|A_S) = Shannon({p_k}) + sum_k p_k h2((1+r_k)/2),  r_k = |(x,y,z)|
#     This is the same mathematics reached by a different route, and it exposes the two pieces
#     the ledger arms never separate:
#         H_CENTRE = Shannon({p_k})      -- THE BOUNDARY'S OWN DATA  (2 bits max)
#         C        = sum_k p_k E_k       -- THE RECORD BEYOND THE BOUNDARY (1 bit max)
#         H_FULL   = H_CENTRE + C
#
# THE IDENTITY BOTH LANES REPORT, RE-DERIVED HERE IN THREE LINES SO IT IS NOT TAKEN ON TRUST:
#   A_env = A_S' (commutant).  For a pure state and A = (+)_k M_{d_k} (x) 1_{m_k}:
#       S(A)  = H(p) + sum p_k E_k        S(A') = H(p) + sum p_k E_k     (block states pure)
#       A v A' = (+)_k M_{d_k} (x) M_{m_k}  =>  S(A v A') = H(p)
#   =>  I(A_S : A_env) = 2 H_FULL - H_CENTRE .   A_Sigma = Z(A_S) is a SUBalgebra of A_S, so
#       I(A_S : A_Sigma) = H_FULL + H_CENTRE - H_FULL = H_CENTRE .
#   =>  Delta_surf = I(A_S:A_env) - I(A_S:A_Sigma) = 2(H_FULL - H_CENTRE) = 2C .
#   The pre-registered PRIMARY FALSIFIER is therefore a function of rho|A_S alone and never
#   reads the environment.  Both lanes found this; it is re-derived independently here.

import numpy as np, itertools, math

# ---------------------------------------------------------------- INHERITED CARRIER
V_N = 8
EDGES = [(0, 7), (0, 1), (0, 2), (1, 2), (1, 3), (2, 3),
         (3, 4), (4, 5), (4, 6), (5, 6), (5, 7), (6, 7)]
L = len(EDGES)
DIMF = 1 << L                                        # 4096

STAR = [0] * V_N
for _i, (_a, _b) in enumerate(EDGES):
    STAR[_a] |= 1 << _i
    STAR[_b] |= 1 << _i

PLAQ = [0b1110, 0b111000, 0b1110000000, 0b111000000000,
        (1 << 0) | (1 << 1) | (1 << 4) | (1 << 6) | (1 << 7) | (1 << 10)]

GRID = [0.05, 0.10, 0.20, 0.30, 0.45, 0.60, 0.80, 1.00, 1.30, 1.70, 2.20, 3.00, 5.00]
G2_SD = math.sqrt(5.0 / 12.0)

S_LINKS = [1, 2, 3]                                  # the region's own links
SIG_LINKS = [0, 4, 5]                                # delta(A) for A = {0,1,2}


def LM(ls):
    m = 0
    for l in ls:
        m |= 1 << l
    return m


S_T = LM(S_LINKS)
SIG_T = LM(SIG_LINKS)


def pop(m):
    return bin(m).count("1")


def dpartial(mask):
    r = 0
    for i in range(L):
        if mask >> i & 1:
            a, b = EDGES[i]
            r ^= (1 << a) ^ (1 << b)
    return r


for _p in PLAQ:
    assert dpartial(_p) == 0, "plaquette is not a cycle"
assert dpartial(S_T) == 0, "S is not a cycle -- W_S undefined"

# ---------------------------------------------------------------- Z-basis config table
_CONF = np.arange(DIMF, dtype=np.int64)
_PARITY = np.zeros(DIMF, dtype=np.int8)
for _b in range(L):
    _PARITY ^= ((_CONF >> _b) & 1).astype(np.int8)


def zsign_full(b):
    """diag of Z^b on the 4096 space, as +-1 float vector."""
    if b == 0:
        return np.ones(DIMF)
    acc = np.zeros(DIMF, dtype=np.int8)
    for j in range(L):
        if b >> j & 1:
            acc ^= ((_CONF >> j) & 1).astype(np.int8)
    return 1.0 - 2.0 * acc.astype(np.float64)


# ---------------------------------------------------------------- physical sector as an isometry
class Sector:
    """The 32-dim physical sector for charge pattern eta, as an explicit 4096 x 32 isometry."""

    def __init__(self, eta):
        eta = list(eta)
        assert len(eta) == V_N and int(np.prod(eta)) == 1
        self.eta = eta
        cut_sign = {}
        for A in range(1 << V_N):
            m, s = 0, 1
            for v in range(V_N):
                if A >> v & 1:
                    m ^= STAR[v]
                    s *= eta[v]
            if m in cut_sign:
                assert cut_sign[m] == s
            cut_sign[m] = s
        assert len(cut_sign) == 128
        self.cut_sign = cut_sign
        seen = np.zeros(DIMF, dtype=bool)
        cols = []
        for c in range(DIMF):
            if seen[c]:
                continue
            col = np.zeros(DIMF)
            for m, s in cut_sign.items():
                col[c ^ m] = s
                seen[c ^ m] = True
            cols.append(col / math.sqrt(128.0))
        self.U = np.array(cols).T                     # 4096 x 32
        self.dim = self.U.shape[1]
        assert self.dim == 32
        assert np.allclose(self.U.T @ self.U, np.eye(32), atol=1e-12)
        self._op = {}

    # ---- operators, restricted to the sector -------------------------------------------------
    def op(self, a=0, b=0):
        """32x32 matrix of the HERMITIAN string O(a,b) = i^{|a&b|} X^a Z^b, restricted."""
        key = (a, b)
        if key in self._op:
            return self._op[key]
        W = self.U * zsign_full(b)[:, None]
        if a:
            W = W[_CONF ^ a, :]
        M = self.U.T @ W
        ph = 1j ** (pop(a & b) & 3)
        M = M.astype(complex) * ph
        assert np.linalg.norm(M - M.conj().T) < 1e-10, ("not hermitian", a, b)
        self._op[key] = M
        return M

    def H(self, g2):
        M = np.zeros((self.dim, self.dim), dtype=complex)
        for p in PLAQ:
            M -= (1.0 / g2) * self.op(0, p)
        for l in range(L):
            M -= g2 * self.op(1 << l, 0)
        return (M + M.conj().T) / 2

    def ground(self, g2):
        ev, evec = np.linalg.eigh(self.H(g2))
        return evec[:, 0].copy(), float(ev[1] - ev[0])

    def eig(self, g2):
        return np.linalg.eigh(self.H(g2))


_SECTOR_CACHE = {}


def sector(charges):
    """charges: iterable of vertices carrying eta = -1.  Must have even size."""
    key = tuple(sorted(charges))
    if key not in _SECTOR_CACHE:
        eta = [-1 if v in key else 1 for v in range(V_N)]
        _SECTOR_CACHE[key] = Sector(eta)
    return _SECTOR_CACHE[key]


def flux_of(charges):
    """X^Sigma = X_0 X_4 X_5 = delta({0,1,2}) = eta_0 eta_1 eta_2 on the physical sector."""
    s = 1
    for v in (0, 1, 2):
        if v in charges:
            s = -s
    return s


# ---------------------------------------------------------------- the record on A_S
def h2(p):
    if p <= 1e-15 or p >= 1 - 1e-15:
        return 0.0
    return float(-p * math.log2(p) - (1 - p) * math.log2(1 - p))


def shannon(ps):
    out = 0.0
    for p in ps:
        if p > 1e-14:
            out -= p * math.log2(p)
    return float(out)


class RecordOps:
    """Concrete operators of the record on S, for one sector."""

    def __init__(self, sec):
        self.sec = sec
        self.c1 = sec.op(LM([1, 2]), 0)               # X_1 X_2  ( = eta_0 X_0 by Gauss )
        self.c2 = sec.op(LM([1, 3]), 0)               # X_1 X_3  ( = eta_1 X_4 by Gauss )
        self.bx = sec.op(1 << 1, 0)                   # X_1
        self.bz = sec.op(0, S_T)                      # W_S = Z_1 Z_2 Z_3
        self.by = 1j * (self.bx @ self.bz)            # Hermitian third direction
        I = np.eye(sec.dim, dtype=complex)
        self.P = {}
        for s1 in (+1, -1):
            for s2 in (+1, -1):
                self.P[(s1, s2)] = ((I + s1 * self.c1) / 2) @ ((I + s2 * self.c2) / 2)
        self.keys = [(+1, +1), (+1, -1), (-1, +1), (-1, -1)]


_REC = {}


def rec_ops(sec):
    k = id(sec)
    if k not in _REC:
        _REC[k] = RecordOps(sec)
    return _REC[k]


def record(sec, psi, tol=1e-12):
    """Full decomposition of rho|A_S.  Returns a dict.

    p[k]      central weights = THE BOUNDARY'S OWN DATA (A_Sigma = Z(A_S), 2 bits)
    bloch[k]  conditional block Bloch vector (x,y,z) inside central sector k
    H_CENTRE  Shannon(p)                       -- boundary data entropy
    C         sum_k p_k E_k                    -- record content BEYOND the boundary
    H_FULL    H_CENTRE + C
    H_BLOCK   vN entropy of the UNCONDITIONED block qubit  (alg{X_1, W_S})
    H_MAG     h2 of <W_S>
    DELTA_SURF = 2C  (the pre-registered primary falsifier)
    """
    R = rec_ops(sec)
    ps, bl, Es = {}, {}, {}
    for k in R.keys:
        v = R.P[k] @ psi
        p = float(np.vdot(v, v).real)
        ps[k] = max(p, 0.0)
        if p <= tol:
            bl[k] = (0.0, 0.0, 0.0)
            Es[k] = 0.0
            continue
        x = float(np.vdot(v, R.bx @ v).real) / p
        y = float(np.vdot(v, R.by @ v).real) / p
        z = float(np.vdot(v, R.bz @ v).real) / p
        bl[k] = (x, y, z)
        r = min(math.sqrt(x * x + y * y + z * z), 1.0)
        Es[k] = h2((1 + r) / 2)
    HC = shannon([ps[k] for k in R.keys])
    C = sum(ps[k] * Es[k] for k in R.keys)
    ux = float(np.vdot(psi, R.bx @ psi).real)
    uy = float(np.vdot(psi, R.by @ psi).real)
    uz = float(np.vdot(psi, R.bz @ psi).real)
    ur = min(math.sqrt(ux * ux + uy * uy + uz * uz), 1.0)
    return {
        "p": ps, "bloch": bl, "E": Es,
        "H_CENTRE": HC, "C": C, "H_FULL": HC + C,
        "H_BLOCK": h2((1 + ur) / 2), "H_MAG": h2((1 + uz) / 2),
        "DELTA_SURF": 2 * C,
        "uncond_bloch": (ux, uy, uz),
    }


def tv(pa, pb, keys):
    return 0.5 * sum(abs(pa[k] - pb[k]) for k in keys)


def dtr_qubit(ba, bb):
    """trace distance between two qubit states given Bloch vectors."""
    d = np.array(ba) - np.array(bb)
    return 0.5 * float(np.linalg.norm(d))


def dtr_record(ra, rb, keys=None):
    """trace distance between rho|A_S of two states, INTRINSIC labelling.
       rho|A_S = (+)_k p_k rho_k  ->  D = 1/2 sum_k || p_k rho_k - q_k sigma_k ||_1 ."""
    if keys is None:
        keys = list(ra["p"].keys())
    tot = 0.0
    for k in keys:
        p, q = ra["p"][k], rb["p"][k]
        A = _qubit(p, ra["bloch"][k])
        B = _qubit(q, rb["bloch"][k])
        w = np.linalg.eigvalsh(A - B)
        tot += 0.5 * float(np.abs(w).sum())
    return tot


def _qubit(p, b):
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    return p * (np.eye(2, dtype=complex) + b[0] * sx + b[1] * sy + b[2] * sz) / 2


# ---------------------------------------------------------------- logging
class Log:
    def __init__(self, path):
        self.path = path
        self.buf = []

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        self.buf.append(s)
        print(s, flush=True)

    def rule(self, t=""):
        self("\n" + "=" * 100)
        self(t)
        self("=" * 100)

    def save(self):
        with open(self.path, "w") as f:
            f.write("\n".join(self.buf) + "\n")
