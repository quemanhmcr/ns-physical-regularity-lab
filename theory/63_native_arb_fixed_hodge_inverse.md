# Native Arb observer for fixed Hodge inverses

The coupled feedback law repeatedly applies two already-validated fixed physical linear operators, `K44^{-1}` and `K66^{-1}`.  Hand-written interval Gaussian elimination is a poor observer for deeply nested evaluations: at 160/256 bits it can inflate small input boxes into enormous enclosures even though the fixed operators are nonsingular and the 512-bit Krawczyk test succeeds with a wide inclusion margin.

Module 135 changes **only the numerical observer for these fixed linear eliminations**.  It uses python-flint's native `arb_mat.solve`, which is designed to solve matrices in ball arithmetic.  The physical bases, Hodge lifts, feedback map, axial symmetry reduction, candidate branch, interval box radii and Krawczyk inclusion criterion are unchanged from module 134.

If the same certificate now survives 160/256/512, the difference between modules 134 and 135 is an observer autopsy: raw Gaussian interval dependency, not a physical obstruction.
