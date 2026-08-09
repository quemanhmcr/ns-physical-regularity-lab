# Impulse-neutral donors and the natural multipole ladder

The previous far-field experiment found that closed vorticity eliminates the monopole and that hydrodynamic impulse is the first surviving remote variable.  The next escape route is to build large local `+I` and `-I` structures so that total impulse remains small.

Do not answer this by imposing a norm on local impulse.  Let the physical cancellation act on the field itself.

For a compact donor with impulse `I`, the leading remote velocity is the impulse dipole

`u_I(x) = [3 e (I.e) - I]/(4 pi r^3)`.

Put a `+I` donor at `+a/2` and a `-I` donor at `-a/2`.  The total impulse is exactly zero.  Expanding their physical fields gives

`u_pair(x) = -(a.grad) u_I(x) + higher terms`.

Thus the leading influence scales as `|I||a|/r^4`.  The translation-invariant first spatial moment of the neutral impulse distribution is the dyad

`M1 = sum c_alpha tensor I_alpha`;

for the pair, `M1 = a tensor I`.  Translation invariance follows because `sum I_alpha = 0`.

For a circular ring, `I = pi Gamma R^2 n`.  A coaxial opposite-circulation pair therefore has exact axial velocity whose first nonzero far term is proportional to `Gamma R^2 a / z^4`.

A second cancellation level is obtained from four closed rings at positions `(-3a/2,-a/2,+a/2,+3a/2)` with signed impulses `(+I,-I,-I,+I)`.  Here both

`sum I_alpha = 0`

and

`sum c_alpha I_alpha = 0`.

The first surviving field is one derivative farther down the hierarchy and scales as `r^-5`.

The physical lesson is narrower than a general theorem but structurally important: cancellation does not hide remote influence for free.  Each cancelled lower moment removes an entire far-field power and promotes the next spatial moment as the variable nature permits a remote observer to see.

The remote tests must therefore distinguish two claims:

1. **Decay ladder:** exact closed-ring stacks exhibit `r^-3 -> r^-4 -> r^-5` when zero, one, then two leading impulse-moment levels are cancelled.
2. **Microscopic-shape forgetting:** an impulse-neutral pair built from different closed loop shapes but with the same `I` and separation `a` converges to the same leading `r^-4` field.

If both survive, the next regularity question is whether a singular cascade can finance an infinite sequence of increasingly high moments when lower moments are forced to cancel by global constraints.
