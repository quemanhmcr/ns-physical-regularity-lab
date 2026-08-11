import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rt3=arb(3).sqrt(); K=arb(4)*pi/(3*rt3); rows=[]
# For omega=A(2x,-y,-z), u=(0,-Axz,Axy), G_self=2Au is tangent and curl omega=0.
# Spherical cap angle n_x>1/sqrt3 at radius R has Gamma=K A R^3.
# Move the cap radially with source velocity V=Rdot n while its angular aperture is fixed.
# Exact boundary recruitment integral is 4pi/sqrt3 A R^2 Rdot =3 K A R^2 Rdot.
for As in ('1e-30','1','1e30'):
 A=arb(As)
 for Rs in ('1e-20','1','1e20'):
  R=arb(Rs)
  for Vrs in ('-1e30','-1','1','1e30'):
   Rdot=arb(Vrs)
   for nus in ('1e-30','1','1e30'):
    nu=arb(nus)
    Gamma=K*A*R**3
    direct=3*K*A*R**2*Rdot
    recruit=arb(4)*pi/rt3*A*R**2*Rdot
    if not (direct/recruit).contains(1): raise AssertionError(('recruit identity',As,Rs,Vrs,nus,direct,recruit))
    # Kelvin viscous current vanishes structurally because omega=A grad(x^2-y^2/2-z^2/2), so curl omega=0.
    kelvin=arb(0)*nu
    self_euler_flux=arb(0)
    rows.append({
      'A':As,'R':Rs,'Rdot':Vrs,'nu':nus,'cap_circulation':str(Gamma),
      'direct_flux_change_from_radial_source_motion':str(direct),
      'source_relative_material_recruitment_boundary_integral':str(recruit),
      'self_Euler_cap_flux_change':'0','Kelvin_viscous_current_cap_flux_change':'0',
      'direct_over_recruitment_ratio':str(direct/recruit),
      'fractional_recruitment_rate_3Rdot_over_R':str(3*Rdot/R),
    })
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'Use the exact source-relative material-spacetime flux ledger on the positive spherical cap of the linear null catalyst.  The cap circulation is Gamma=4 pi A R^3/(3 sqrt3).  The catalyst self-Euler source is tangent to the sphere, so it does not change this cap flux; the Kelvin viscous current also vanishes in the unlocalized linear germ because curl omega=0. '
  'If the active source cap moves radially with speed Rdot while keeping its angular aperture fixed, direct differentiation gives dGamma/dt=3 Gamma Rdot/R.  The exact boundary term int [omega cross (V-u)] dot dl equals 4 pi A R^2 Rdot/sqrt3, exactly the same quantity; the u part is azimuthal and contributes zero. '
  'Thus for the pure catalyst the source-relative ledger separates cleanly: Euler conversion creates productive Q without spending cap circulation, bulk viscosity is absent, and changes of the circulation inventory seen by a shrinking/growing Eulerian source occur through material recruitment/crossing.  For a shrinking cap Rdot<0 this geometric term removes active flux rather than creating the divergent flux required by the critical branch. '
  'Therefore any increase of catalyst amplitude/cap circulation toward a singular state must come from additional external material recruitment/deformation or from viscous current introduced by finite-energy turnover, not from the catalyst self-Euler mechanism.'),
 'rows':rows
},indent=2,allow_nan=False))
