"""rf2ls: the linear system of a partial fraction decomposition, built column by column.

Dependency-free Python (standard library only) companion to the MATLAB
routine RF2LS.m. All arithmetic is exact (fractions.Fraction).

Polynomials are lists of coefficients, highest power first:
    s^2 + 4  ->  [1, 0, 4]
The denominator is given already factored, as a list of
(factor, multiplicity) pairs, where each factor is monic and linear or an
irreducible quadratic:
    (s + 1) (s^2 + 4)    ->  [([1, 1], 1), ([1, 0, 4], 1)]
    s (s + 1)^2          ->  [([1, 0], 1), ([1, 1], 2)]

Example
-------
>>> M, n, labels, lam = rf2ls([3, 1], [([1, 1], 1), ([1, 0, 4], 1)])
>>> [[str(x) for x in row] for row in M]
[['1', '1', '0'], ['0', '1', '1'], ['4', '0', '1']]
>>> [str(x) for x in lam]
['-2/5', '2/5', '13/5']

Author: Carlos Alberto Pereyda-Pierre. License: MIT.
"""

from fractions import Fraction

__all__ = ["rf2ls", "polyvec", "format_decomposition"]


def _poly(p):
    return [Fraction(c) for c in p]


def _mul(p, q):
    r = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i + j] += a * b
    return r


def _prod(polys):
    r = [Fraction(1)]
    for p in polys:
        r = _mul(r, p)
    return r


def polyvec(p, delta):
    """Coefficient vector of p in the basis {s^(delta-1), ..., s, 1}."""
    p = _poly(p)
    while len(p) > 1 and p[0] == 0:          # drop leading zeros
        p = p[1:]
    if len(p) > delta:
        raise ValueError(f"degree {len(p) - 1} exceeds delta - 1 = {delta - 1}")
    return [Fraction(0)] * (delta - len(p)) + p


def _solve(M, b):
    """Gauss-Jordan elimination with exact fractions."""
    n = len(b)
    A = [list(row) + [b[i]] for i, row in enumerate(M)]
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            raise ValueError("singular matrix: check that the factors are distinct and irreducible")
        A[c], A[piv] = A[piv], A[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c] / A[c][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [A[i][n] / A[i][i] for i in range(n)]


def rf2ls(N, factors):
    """Build and solve M * lambda = n for the proper rational function N/D.

    Parameters
    ----------
    N : list of coefficients of the numerator (highest power first).
    factors : list of (factor, multiplicity) giving D = prod factor**multiplicity.

    Returns
    -------
    M : list of rows of the delta-by-delta matrix (one column per unknown).
    n : coefficient vector of N.
    labels : description of the unknown attached to each column.
    lam : the partial fraction coefficients, in the order of `labels`.
    """
    for f, m in factors:
        if len(f) not in (2, 3) or Fraction(f[0]) != 1 or m < 1:
            raise ValueError("each factor must be monic, of degree 1 or 2, with multiplicity >= 1")
    delta = sum((len(f) - 1) * m for f, m in factors)
    Nr = polyvec(N, len(N))
    while len(Nr) > 1 and Nr[0] == 0:
        Nr = Nr[1:]
    if len(Nr) - 1 >= delta:
        raise ValueError("the rational function must be proper: deg N < deg D")

    cols, labels = [], []
    for idx, (f, m) in enumerate(factors):
        f = _poly(f)
        R = _prod([_poly(g) for j, (g, mm) in enumerate(factors) if j != idx for _ in range(mm)])
        for k in range(1, m + 1):
            Q = _mul(R, _prod([f] * (m - k)))            # D / f^k
            if len(f) == 2:
                cols.append(Q)
                labels.append(("A", f, k))
            else:
                cols.append(_mul([Fraction(1), Fraction(0)], Q))   # s * D / f^k
                labels.append(("B", f, k))
                cols.append(Q)
                labels.append(("C", f, k))

    vecs = [polyvec(c, delta) for c in cols]
    M = [[v[i] for v in vecs] for i in range(delta)]
    n = polyvec(N, delta)
    lam = _solve(M, n)
    return M, n, labels, lam


def _fmt_poly(p):
    p = _poly(p)
    deg = len(p) - 1
    terms = []
    for i, c in enumerate(p):
        e = deg - i
        if c == 0:
            continue
        mono = "" if e == 0 else ("s" if e == 1 else f"s^{e}")
        coef = str(c) if (mono == "" or abs(c) != 1) else ("-" if c < 0 else "")
        terms.append(f"{coef}{'*' if mono and coef not in ('', '-') else ''}{mono}")
    return " + ".join(terms).replace("+ -", "- ") or "0"


def _c(x):
    return str(x) if (x.denominator == 1 and x >= 0) else f"({x})"


def format_decomposition(labels, lam):
    """Human-readable sum of partial fractions."""
    parts, i = [], 0
    while i < len(labels):
        kind, f, k = labels[i]
        d = f"({_fmt_poly(f)})" + (f"^{k}" if k > 1 else "")
        if kind == "A":
            parts.append(f"{_c(lam[i])}/{d}")
            i += 1
        else:
            parts.append(f"({_c(lam[i])}*s + {_c(lam[i + 1])})/{d}")
            i += 2
    return " + ".join(parts)


if __name__ == "__main__":
    examples = [
        ("(5s+7)/((s+1)(s+2))", [5, 7], [([1, 1], 1), ([1, 2], 1)]),
        ("(4s+7)/(s(s+1)^2)", [4, 7], [([1, 0], 1), ([1, 1], 2)]),
        ("(3s+1)/((s+1)(s^2+4))", [3, 1], [([1, 1], 1), ([1, 0, 4], 1)]),
        ("(22s+33)/((s+98)^2(s+76))", [22, 33], [([1, 98], 2), ([1, 76], 1)]),
    ]
    for name, N, factors in examples:
        M, n, labels, lam = rf2ls(N, factors)
        print(f"\nF(s) = {name}")
        print("M =", [[str(x) for x in row] for row in M])
        print("n =", [str(x) for x in n])
        print("F(s) =", format_decomposition(labels, lam))
