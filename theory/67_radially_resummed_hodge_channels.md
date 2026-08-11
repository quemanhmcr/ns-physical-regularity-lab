# Radially resummed Hodge channels

## THINK — polynomial degree was partly an observer coordinate

A toroidal angular vorticity channel should not be represented as unrelated copies `r^q T_l`.  Let

`omega = a(r) x cross grad H_l`,

where `H_l` is a solid harmonic of angular degree `l` and `a(r)` is an arbitrary smooth radial profile.

Seek its tangent Hodge velocity as

`u = A(r) grad H_l + B(r) H_l x`.

The exact div-curl equations give

`A'/r-B=a`,

and

`A''+(2l+2)A'/r = r a' +(l+3)a`.

Define

`I(r)=integral_0^r s^(2l+2) a(s) ds`.

Regularity at the center gives

`A'(r)=r a(r)-l r^(-2l-2) I(r)`,

`B(r)=-l r^(-2l-3) I(r)`.

Tangency at the physical source sphere `r=L` fixes the remaining harmonic constant.  The degree-`l-1` harmonic companion is

`C_l[a] grad H_l`,

where

`C_l[a]=-(l+1)/(2l+1) integral_0^L r a(r)[1-(r/L)^(2l+1)] dr`.

This is the native radial Hodge screen for angular channel `l`.

## The original kernel is the l=2 member

For `l=2`,

`C_2[a]=-(3/5) integral_0^L r a(r)[1-(r/L)^5] dr`.

Thus the kernel `1-(r/L)^5` that first appeared in the strain Hodge transaction microscope is not an isolated exponent.  It is the quadrupole member of the general family

`1-(r/L)^(2l+1)`.

The same physical boundary problem selects the screen at every angular order.

## Consequence for the servo hierarchy

Every radial Taylor copy of the same toroidal `T_l` channel shares the same lower harmonic direction `grad H_l`.  Its entire lower-level influence is summarized by the single screened moment `C_l[a]`.

Therefore the polynomial triangular matrices are coordinates of a simpler object:

`angular Hodge channel + radial profile -> one screened harmonic feedback amplitude`.

This is the natural language in which further coupled maintenance should be formulated.  Polynomial homogeneity remains useful for exact computation, but it is demoted as a fundamental physical hierarchy.
