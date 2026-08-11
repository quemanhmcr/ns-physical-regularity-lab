# Screened radial feedback routing

For a toroidal angular channel

`omega_l = a_l(r) x cross grad H_l`,

the exact tangent Hodge lift contains one lower harmonic velocity component

`U_low = C_l[a_l] grad H_l`,

with

`C_l[a]=-(l+1)/(2l+1) integral_0^L r a(r)[1-(r/L)^(2l+1)] dr`.

Therefore its action on the base productive vorticity is exactly

`[omega_2,U_low] = C_l[a] [omega_2,grad H_l]`.

The lower response has homogeneous degree `l`, independent of the radial Taylor order used to represent `a(r)`.

Thus the physical routing diagram is not

`r^0 T_l -> one lower response`, `r^2 T_l -> another`, ...

but

`T_l[a(r)] -> screened scalar C_l[a] -> one angular lower-response channel`.

The radial kernel `ker C_l` is infinite-dimensional.  A simple explicit element at `L=1` is

`a_l^0(r)=1-[2(2l+5)/(2l+3)] r^2`,

for which `C_l[a_l^0]=0` although the vorticity field is nonzero.

This forces an important reinterpretation of the coupled degree-four/six Krawczyk result.  The certified 3/5/9-dimensional fixed point is a locally unique equilibrium of the **screened feedback moments inside the calibrated polynomial response manifold**.  It is not a uniqueness theorem for all radial Hodge profiles.  Internal `ker C_l` radial shapes remain available to Navier-Stokes and may modify higher nonlinear emission without altering lower harmonic feedback.

This radial-silent freedom must be attacked before any claim that higher emission creates an unavoidable infinite servo cost.
