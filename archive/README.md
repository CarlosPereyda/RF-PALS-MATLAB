# Archived: PALS2RF (do not use)

`PALS2RF.m` tried to go in the opposite direction: from an arbitrary square
linear system back to a rational function, by reading each column as a
polynomial and taking their least common multiple as the denominator.

That reverse construction is **not correct in general**, so it is kept here
only for the record and is no longer part of the toolbox.

Why it fails:

* **Quadratic factors.** The column of a coefficient of `s` is `s*D/f^k`,
  so the extra factor `s` enters the least common multiple. For
  `(3s+1)/((s+1)(s^2+4))` the columns are `s^2+4`, `s(s+1)` and `s+1`, whose
  lcm is `s(s+1)(s^2+4)`: a denominator of the wrong degree.
* **A single repeated linear factor.** For `1/(s+1)^2` the columns are
  `s+1` and `1`, whose lcm is `s+1`, so the reconstructed function is
  `1/(s+1)`.
* **Arbitrary systems.** The identity `F = sum(lambda_r / Q_r)` with
  `Q_r = D/P_r` holds for *any* invertible matrix, but the `Q_r` are in
  general not powers of irreducible factors, so the result is not a partial
  fraction decomposition (the 3-by-3 identity matrix gives `D = s^2` and an
  improper function).
* **Lost information.** A matrix alone does not say which columns belong
  to `B` (coefficient of `s`) and which to `C` (constant) unknowns.

The forward direction, implemented in `src/RF2LS.m`, is always valid: every
proper rational function corresponds to exactly one square system, and that
system has exactly one solution.
