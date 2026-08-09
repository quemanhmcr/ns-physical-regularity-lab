# Harmonic-vortical split: a source-free law without a source model

The remote-donor program should not rely on a literal vorticity-free cavity in viscous Navier-Stokes, because viscosity destroys compact support. The more natural local object is the Hodge split of the actual velocity inside a ball `B_d(x0)`.

Let `h = grad phi` solve

`Delta phi = 0` in `B_d`,

with `partial_n phi = u.n` on the boundary (and an arbitrary additive constant fixed by zero mean). Define `v=u-h`. Then

`div h = div v = 0`, `curl h = 0`, and `v.n=0` on the boundary.

The split is kinetic-energy orthogonal:

`integral_B h.v = 0`,

so

`E_B(u) = E_B(h)+E_B(v)`.

Moreover `curl v = curl u = omega` in the ball: all local vorticity resides in the vortical/tangential component, while `h` is the irrotational potential-deformation component selected by the actual boundary flux.

Now expand the harmonic potential around the center. The degree-two harmonic term is

`phi_2(x)=1/2 x.S_h(0)x`,

where `S_h(0)` is symmetric and traceless. Different homogeneous harmonic degrees are orthogonal in the Dirichlet inner product on a ball. Therefore every translation or higher harmonic mode only adds kinetic energy beyond the degree-two contribution. Since

`integral_{B_d} x_i x_j dx = (4pi d^5/15) delta_ij`,

we obtain the exact floor

`E_h(B_d) >= (2pi/15) |S_h(x0)|_F^2 d^5`.

Equality is attained by pure linear strain `h(x)=S_h(x0)(x-x0)`.

This is a sharper and more natural version of the donor-horizon idea. It does not ask where a remote vortex donor sits or what shape it has. It asks how much of the local strain is carried by the irrotational harmonic component versus the vortical component.

If the harmonic part carries strain `s_h`, finite total energy `E0` implies

`d <= [15 E0/(2pi s_h^2)]^(1/5)`.

If a large total strain is not carried by `h`, then it must reside in `v`, the component that contains the local vorticity and is therefore the correct object for the vorticity-direction/viscous-geometry analysis.

This gives an exact physical dichotomy:

- potential/harmonic deformation pays a `d^5` kinetic-energy floor;
- vortical deformation cannot be outsourced to a source model and must be analyzed through local vorticity geometry.

The next task is to determine whether the vortical component's strain at the center can remain singularly efficient while avoiding the exact directional-viscous term `-nu |omega| |grad xi|^2`.
