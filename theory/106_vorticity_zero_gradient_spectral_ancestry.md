# Spectral ancestry law at a vorticity zero

Let

`A=grad u`, `B=grad omega`.

Differentiate the vorticity equation along a material trajectory. At a point where `omega=0`,

`D_t B=[A,B]+nu grad Delta omega`.

The non-affine Euler term `(grad A)omega` vanishes exactly.

Therefore

`D_t tr(B^k)=k nu tr(B^(k-1) grad Delta omega)`.

In Euler all spectral invariants of the linear vorticity gradient are frozen. Strain may conjugate `B`, producing large coordinate entries in a nonnormal frame, but it cannot change its eigenvalue invariants.

For the symmetric linear null catalyst this means its actual spectral strength cannot be increased by the same reusable Euler mechanism that lets it emit productive strain. The only channels are:

1. viscous third-derivative current;
2. switching/recruitment to another material vorticity zero carrying a different `B` spectrum.

This is a much sharper ancestry statement than a norm estimate: it identifies which part of the local Taylor germ is Euler-reusable and which part can only mutate viscously.
