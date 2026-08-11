# Simple-vorticity-zero genealogy is viscous

Let `z(t)` be a simple vorticity zero:

`omega(z(t),t)=0`, `B=grad omega`, `det B !=0`.

Differentiating the zero condition with the exact Navier-Stokes vorticity equation gives

`B(V_z-u)+nu Delta omega=0`.

Hence

`V_z-u=-nu B^-1 Delta omega`.

A simple vorticity zero is therefore material in Euler. Its relative drift through the fluid is purely viscous.

Following the zero worldline,

`dB/dt=[A,B]+nu grad Delta omega +(V_z-u).grad B`.

Substituting the drift law makes the last two terms explicitly proportional to viscosity. Consequently all spectral invariants of `B` change only by viscous zero-genealogy current; the Euler commutator changes eigenvectors but not eigenvalues.

Thus recruiting a *fresh simple zero* does not evade the ancestry law. The only new geometric possibility is a degenerate zero, `det B=0`, where isolated zero branches merge into a zero curve/surface and the inverse-B drift law loses transversality.
