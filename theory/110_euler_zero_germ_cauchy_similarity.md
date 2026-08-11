# Euler Cauchy similarity of a vorticity-zero germ

Euler's Cauchy formula is

`omega(t,X(a,t))=F(a,t) omega0(a)`, `F=grad_a X`, `det F=1`.

At an initial vorticity zero `omega0(a)=0`, differentiate with respect to the current coordinate. The term involving `grad_a F` is multiplied by `omega0` and vanishes. Therefore

`B(t)=F B0 F^-1`, `B=grad omega`.

So at every Euler zero ancestry:

- eigenvalues of `B` are fixed;
- `tr(B^k)` and `det B` are fixed;
- rank and Jordan structure are fixed;
- for a degenerate zero manifold, `ker B` is transported by `F`.

Degeneracy can be deformed but not upgraded spectrally by Euler. A source that samples stronger and stronger zero germs must either access stronger germs already present in the initial zero set or rely on viscosity to mutate/create them.
