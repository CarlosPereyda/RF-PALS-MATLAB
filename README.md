# RF-PALS-MATLAB

MATLAB implementation of the RF--PALS correspondence framework introduced in:

"A Structured Correspondence Between Rational Functions and Polynomially Associated Linear Systems"

This repository contains symbolic and matrix-oriented MATLAB routines for constructing Polynomially Associated Linear Systems (PALS) from proper rational functions and reconstructing associated rational functions from linear systems.

---

## Repository Structure

RF-PALS-MATLAB/
│
├── src/
│   ├── RF2PALS.m
│   ├── PALS2RF.m
│   └── polyvec.m
│
├── examples/
│
└── README.md

---

## Main Functions

### RF2PALS.m

Constructs the Polynomially Associated Linear System (PALS) induced by a proper rational function.

### PALS2RF.m

Reconstructs an associated rational function from a Polynomially Associated Linear System.

### polyvec.m

Converts symbolic polynomials into coefficient vectors in the ordered monomial basis.

---

## Requirements

- MATLAB
- Symbolic Math Toolbox

Tested using MATLAB symbolic capabilities for:
- polynomial factorization,
- symbolic simplification,
- least common multiple operations,
- and matrix-based reconstruction procedures.

---

## Conceptual Overview

The RF--PALS framework establishes a structured correspondence between:

- proper rational functions,
- generated polynomials,
- coefficient-vector representations,
- and Polynomially Associated Linear Systems (PALS).

Within this framework:
- rational functions induce structured linear systems,
- and consistent linear systems may reconstruct associated rational functions.

The implementation emphasizes:
- generated polynomial structure,
- coefficient-vector association,
- and matrix-oriented representation,
rather than traditional procedural partial fraction decomposition.

---

## Author

Carlos Alberto Pereyda-Pierre

---

## License

This project is distributed under the MIT License.