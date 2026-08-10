# Capacity, transaction Reynolds number, and material-spacetime ancestry

## THINK — do not force a deformed collar back into a radius

The straight annulus revealed a harmonic circulation mode

`E_h = Gamma^2 ell log(b/a)/(4 pi)`.

The factor `log(b/a)` is not fundamentally a radius formula.  It is the reciprocal of the physical conductance of the vorticity-free collar cross-section.

Let `A` be any smooth doubly connected planar collar cross-section.  Let `phi` solve

`Delta phi = 0` in `A`,

`phi=0` on the inner boundary, `phi=1` on the outer boundary.

Define its Dirichlet capacity

`Cap(A) = integral_A |grad phi|^2 dA`.

Let `J` rotate planar vectors by ninety degrees.  Then

`h_1 = J grad(phi) / Cap(A)`

is divergence-free, curl-free, tangent to both collar boundaries, and has unit circulation around the noncontractible cycle.  Its squared L2 norm is

`integral_A |h_1|^2 dA = 1/Cap(A)`.

Therefore the unique harmonic circulation mode with circulation `Gamma` has exact kinetic energy per unit vortex-line length

`e_h(A,Gamma) = Gamma^2/[2 Cap(A)]`.

This is a shape-independent Hodge/Dirichlet duality in every planar doubly connected cross-section.  The cylinder formula is only the circular representative.

For a circular annulus,

`Cap = 2 pi/log(b/a)`.

For a flat periodic strip of circulation length `P` and breach thickness `d`,

`Cap = P/d`.

For a confocal elliptical annulus with conformal coordinates `mu1<mu<mu2`,

`Cap = 2 pi/(mu2-mu1)`.

In all three cases

`2 e_h Cap = Gamma^2`.

The natural static collar observable is therefore not Euclidean thickness alone but the **circulation conductance/capacity of the actual collar geometry**.

## PREDICT — cheap occupancy and easy penetration are dual, but gap is not the whole story

If a persistent ancestry keeps `Gamma` fixed while its harmonic occupancy density becomes small, the exact identity forces `Cap(A)` to become large.

Large capacity means a unit potential difference can drive a large harmonic flux from the inner to outer collar boundary.  It is therefore a geometric measure of ease of transverse communication.

But capacity is not the same as minimum Euclidean gap.

In a flat periodic strip,

`Cap=P/d`.

Keeping `d` fixed while sending the circulation-loop extent `P->infinity` makes `Cap->infinity` and `e_h->0` without narrowing the gap.  The escape spends circumferential/interface extent instead.

Thus the correct static statement is

> cheap circulation occupancy forces large transverse conductance, not necessarily a small pointwise gap.

A future contradiction must decide whether large conductance is realized by near-contact, large interface extent, folding/corrugation, or another actual geometry.  Replacing capacity by a single distance would throw away precisely the structure nature is using.

## ATTACK I — the transaction Reynolds number was already hiding in the Hodge microscope

For the signed productive shell transaction

`q_e(r)=e.Q(r)e`,

define its circulation-dimensional form

`Gamma_Q,e(r)=r^2 q_e(r)`.

The only viscosity-normalized dimensionless number built from this transaction at the same physical radius is

`R_Q,e(r) = Gamma_Q,e(r)/nu = q_e(r) r^2/nu`.

After taking absolute value only for the clock interpretation,

`|R_Q,e| = [r^2/nu] / [1/|q_e|]`.

Therefore the productive circulation Reynolds number is **exactly** the ratio

`viscous crossing time / productive strain time`.

There are not two independent ledgers here.  The circulation currency and the local stretch-versus-diffusion clock are the same dimensionless transaction.

The regimes are physically distinct:

- `|R_Q,e| << 1`: the shell diffusion clock is shorter than one productive e-folding time;
- `|R_Q,e| ~ 1`: productive and viscous clocks meet;
- `|R_Q,e| >> 1`: the signed transaction can act faster than diffusion across that radius, but it carries a viscosity-scale-or-larger productive circulation transaction.

This is an identity, not a regularity theorem.  A transient or continually regenerated donor may still act when the local clock is small.  But such regeneration is then part of the ancestry-renewal problem rather than a free persistence branch.

### Finite-filament calibration

For the exact finite straight filament test,

`E(alpha)=4 pi d^2 s/Gamma`,

where `s` is its target shear strain and `alpha=L/(2d)`.

Thus

`chi = s d^2/nu = E(alpha)/(4 pi) * (Gamma/nu)`.

Again the strain/diffusion clock and circulation Reynolds number differ only by the actual geometric efficiency.  An order-one efficient donor with `Gamma/nu << 1` necessarily has `chi << 1`: its own scale diffuses faster than the strain time unless geometry continually rebuilds it.

This gives a representation-level explanation for the exact Oseen survival gate rather than treating that gate as an isolated canonical fact.

## ATTACK II — planar capacity is exact where nature supplies slices, but 3D need not supply them

It would be tempting to foliate every vortex collar by planar cross-sections and apply the capacity duality slice by slice.  That is not geometrically legitimate in general.

Let `xi` be the unit vorticity direction and let

`alpha = xi . dx`.

The planes orthogonal to `xi` are locally integrable into surfaces only when the Frobenius condition holds:

`alpha wedge d alpha = 0`,

or equivalently

`xi . curl(xi) = 0`.

The exact helical Beltrami calibration has

`xi=(cos(kz),sin(kz),0)`,

`curl xi = -k xi`,

so

`xi.curl(xi)=-k != 0`.

There is no family of surfaces everywhere orthogonal to the vorticity direction.  Any proof that silently assumes such cross-sections has excluded genuine 3D helical geometry before the physics begins.

Therefore:

- use the exact 2D capacity duality when a real collar-slice structure exists;
- use the full 3D harmonic cohomology class when it does not;
- do not manufacture cross-sections merely to recover a preferred formula.

## ATTACK III — the true ancestry object is a conserved flux in material spacetime

The static collar still does not solve renewal and double counting.  Navier-Stokes itself gives a more intrinsic spacetime structure.

Let `Omega=omega.dA` denote the vorticity flux 2-form and let `S(t)` be any material surface.  The exact vorticity equation gives

`d/dt integral_{S(t)} Omega`

` = nu integral_{S(t)} Delta omega . n dA`

` = -nu integral_{partial S(t)} curl(omega).dl`.

Define the viscous vorticity-flux current 1-form

`j = nu curl(omega).dl`.

Pull both objects back by the material flow `X_t`:

`Omega_hat(t)=X_t^* Omega(t)`,

`j_hat(t)=X_t^* j(t)`.

Then the exact ancestry law is

`partial_t Omega_hat + d j_hat = 0`.

Equivalently, on material label-space times time, the spacetime 2-form

`F_ancestry = Omega_hat - dt wedge j_hat`

is closed:

`d_4 F_ancestry = 0`.

This is a cleaner statement than saying that viscosity "destroys" vortex lines.  Vorticity flux is redistributed through a closed spacetime current.  A material ancestry can persist vertically in time or flow sideways through a viscous current, but it cannot disappear by observer relabeling.

For any material spacetime control volume, internal current interfaces cancel by Stokes.  This gives a natural non-double-counting ledger: only flux through the external spacetime boundary changes the ancestry inventory.

## ATTACK IV — exact Oseen spacetime telescope

Nested fixed-radius circles in the Lamb-Oseen vortex are material.  Let

`G_R(t)=Gamma[1-exp(-R^2/(4 nu t))]`

be the ancestry inventory inside radius `R`, and define the outward viscous circulation current

`J_R(t)=-dG_R/dt`.

For a material annulus `R_i<R<R_{i+1}`,

`d/dt [G_{R_{i+1}}-G_{R_i}] = J_{R_i}-J_{R_{i+1}}`.

Integrating over a time slab gives

`Delta inventory_i = Loss(R_i)-Loss(R_{i+1})`.

Summing neighboring annuli cancels every internal current exactly.  The total ancestry change depends only on the outer spacetime boundary.

This is the bookkeeping property that the earlier Eulerian donor-switch language lacked.  It is not yet a cost inequality, but it identifies the correct object on which any nonreusable cost must live.

## AUTOPSY — what survives the deformation attacks?

### Radius-based collar inductance

Killed as a universal representation.  The exact invariant in a planar collar is capacity, and different shapes realize the same or larger capacity in physically different ways.

### Minimum-gap-only breach law

Killed.  A periodic strip can make capacity arbitrarily large by increasing circulation-loop/interface extent while keeping the transverse gap fixed.

### Universal planar slicing of a 3D vortex tube

Killed by Frobenius.  Helical Beltrami geometry has `xi.curl xi != 0`, so orthogonal cross-sections do not exist.

### Static capacity as the final renewal cost

Not promoted.  Capacity measures spatial conductance, but the single-loop attack already showed that renewal is parabolic.  The finite transfer must be followed in material spacetime.

### Material-spacetime ancestry conservation

Promoted as an exact Navier-Stokes structure.  It is independent of instantaneous vortex-line labeling and telescopes across internal control surfaces.

### Small productive circulation

Sharpened rather than killed universally.  The Hodge transaction itself has the exact clock `|Gamma_Q|/nu`.  If that number is small, persistence at the same scale loses the race to diffusion; any successful singular mechanism must then use rapid regeneration/renewal rather than pretending the donor is long-lived.

## PROMOTE / KILL

### PROMOTE

1. **Planar Hodge-capacity duality:** for every smooth doubly connected planar collar, `e_h=Gamma^2/(2 Cap)` exactly.
2. **Capacity, not radius, is the natural static transverse conductance observable.**
3. **Productive transaction Reynolds number:** `R_Q=q_e r^2/nu=Gamma_Q,e/nu` is exactly the local viscous-clock/productive-time ratio.
4. **Frobenius gate:** a 2D cross-sectional reduction is physically admissible only where the vorticity-normal plane field is integrable; helical geometry proves this is not automatic.
5. **Material-spacetime ancestry conservation:** `partial_t Omega_hat + d j_hat=0`, equivalently `d_4 F_ancestry=0`.
6. **Spacetime telescoping is the natural anti-double-counting structure.**

### KILL / DO NOT PROMOTE

1. Do not replace a deformed collar by an effective radius before identifying its actual capacity/cohomology.
2. Do not infer small Euclidean breach distance from large capacity alone; interface extent can also raise capacity.
3. Do not impose orthogonal vortex cross-sections in genuine 3D geometry when Frobenius fails.
4. Do not treat instantaneous viscous current as a finite renewal cost; thin concentration defeats that observer.
5. No universal parabolic cost for arbitrary 3D collar breach is proved yet.
6. No regularity contradiction is claimed.

## Next frontier — productive ancestry worldsheets

The next natural object is no longer a sequence of Eulerian donors.  It is a material-spacetime worldsheet carrying the circulation flux that actually supplies signed Hodge transactions.

For each productive event, ask:

- which portion of the closed ancestry 2-form `F_ancestry` carries the relevant signed transaction `Gamma_Q,e`;
- whether that portion persists to the next amplification event or crosses a viscous spacetime boundary;
- what spatial capacity/cohomology the surrounding collar has while it persists;
- what spacetime extent is required when the transaction Reynolds number drops to order one or below.

The target is a genuine persistence-or-renewal theorem in material spacetime.  Persistence should be obstructed by finite instantaneous occupancy and frozen-flux thinning; renewal should be measured by flux crossing of a closed spacetime ancestry form, so internal relabeling and double counting disappear automatically.
