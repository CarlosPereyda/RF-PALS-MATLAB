# RF2LS: partial fractions, column by column

MATLAB and Python routines that build the linear system of a partial fraction
decomposition **column by column**, as a change of basis.

For a proper rational function `F(s) = N(s)/D(s)`, each unknown coefficient
of the decomposition has an *associated polynomial*: `D(s)` divided by the
denominator of its fraction, multiplied by `s` when the unknown is the
coefficient of `s` in a quadratic numerator `B s + C`. The coefficient vectors
of these polynomials, in the monomial basis `{s^(δ-1), ..., s, 1}` with
`δ = deg D`, are the columns of a square matrix `M`. Then

```
M * lambda = n,      n = coefficients of N(s),
```

where `lambda` holds the partial fraction coefficients. The associated
polynomials form a basis of the polynomials of degree less than `δ`, so `M` is
always invertible: every proper rational function corresponds to exactly one
such system, and that system has exactly one solution. `M` is the
change-of-basis matrix from the partial fraction basis to the monomial basis.

The system is the same one obtained by comparing coefficients, but it is
written one column at a time, expanding only short numerical polynomials
instead of the full symbolic identity.

## Example

```
F(s) = (3s + 1) / ((s + 1)(s^2 + 4)) = A/(s + 1) + (B s + C)/(s^2 + 4)

associated polynomials:  A -> s^2 + 4,   B -> s(s + 1),   C -> s + 1

[1 1 0] [A]   [0]
[0 1 1] [B] = [3]        A = -2/5,  B = 2/5,  C = 13/5
[4 0 1] [C]   [1]
```

## Contents

```text
src/
  RF2LS.m          builds and solves M*lambda = n from N(s) and D(s)
  polyvec.m        coefficient vector of a polynomial in the monomial basis
examples/
  example_RF2LS.m  the examples of the article, including a symbolic exam item
python/
  rf2ls.py         dependency-free Python version (exact fractions)
  test_rf2ls.py    tests (python -m unittest)
archive/
  PALS2RF.m        earlier reverse construction, kept for the record only;
                   see archive/README.md for why it is not valid in general
```

## MATLAB

Requires the Symbolic Math Toolbox.

```matlab
addpath src
syms s
[M,n,P,den,numPF,lambda,Fpf] = RF2LS(3*s + 1, (s + 1)*(s^2 + 4), s);
disp(M), disp(lambda), pretty(Fpf)
```

Coefficients may be symbolic parameters, which is useful for preparing
personalized exam items:

```matlab
syms a b c d
[M,n,~,~,~,lambda] = RF2LS(a*s + b, (s + c)^2*(s + d), s);
simplify(lambda)   % equivalent to [(a*d - b)/(c - d)^2; (b - a*c)/(d - c); (b - a*d)/(c - d)^2]
```

`help RF2LS` describes all inputs, outputs and the ordering of the unknowns.

## Python

Standard library only. The denominator is given already factored:

```python
from rf2ls import rf2ls, format_decomposition
# (3s + 1) / ((s + 1)(s^2 + 4))
M, n, labels, lam = rf2ls([3, 1], [([1, 1], 1), ([1, 0, 4], 1)])
print(format_decomposition(labels, lam))
# (-2/5)/(s + 1) + ((2/5)*s + (13/5))/(s^2 + 4)
```

Run `python rf2ls.py` for the examples and `python -m unittest` for the tests
(from the `python/` folder).

## Scope and limitations

* Proper rational functions with real coefficients; the denominator must
  factor into linear and quadratic factors. In MATLAB, `factor` works over the
  rationals: a quadratic such as `s^2 - 2` is kept as a quadratic block, which
  still gives a valid (not fully split) decomposition.
* The routines are meant for teaching and for the moderate degrees found in
  coursework. They are not optimized numerical algorithms.

## Citation

If you use these routines, please cite the accompanying article (manuscript in
preparation) and this repository; see `CITATION.cff`.

## License

MIT. Author: Carlos Alberto Pereyda-Pierre.
