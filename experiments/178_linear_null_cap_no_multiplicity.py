import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
# For omega=A(2x,-y,-z), a vortex line in x>0 has
# x'=2Ax, rho'=-A rho, hence C=x rho^2 is invariant and x is strictly monotone.
# Along one line r^2=x^2+C/x.  The positive spherical cap condition omega.n>0 is
# 2x^3>C, exactly the branch where dr^2/dx=2x-C/x^2>0.
for Cs in ('1e-60','1','1e60'):
 C=arb(Cs)
 xturn=(C/2)**(arb(1)/3)
 for ms in ('1.000001','2','1e6'):
  m=arb(ms);x=m*xturn;rho2=C/x;r2=x*x+rho2
  dr2dx=2*x-C/(x*x)
  radial_sign=2*x*x-rho2
  cap_margin=3*x*x-r2
  if not (dr2dx.lower()>0 and radial_sign.lower()>0 and cap_margin.lower()>0):
   raise AssertionError(('outward branch monotonicity',Cs,ms,x,dr2dx,radial_sign,cap_margin))
  rows.append({'invariant_C_xrho2':Cs,'x_over_turning_x':ms,'turning_x':str(xturn),'x':str(x),'rho_squared':str(rho2),'r_squared':str(r2),'dr2_dx':str(dr2dx),'outward_radial_flux_margin_2x2_minus_rho2':str(radial_sign),'positive_cap_margin_3x2_minus_r2':str(cap_margin)})
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'The exact axisymmetric linear null catalyst has vortex-line equations x_prime=2A x, rho_prime=-A rho.  In x>0 the oriented coordinate x is strictly monotone and C=x rho^2 is a first integral.  Hence each line is described by r^2(x)=x^2+C/x, with a single radial turning point 2x^3=C. '
  'The positive outward cap of any centered sphere is exactly the branch omega.n>0, equivalently 2x^3>C.  On that branch dr^2/dx>0 strictly.  Therefore a given oriented vortex line can intersect the positive cap of a given centered sphere at most once. '
  'The positive cap circulation of the pure linear catalyst is therefore not a winding multiplicity count of one lineage at a fixed time: it is the aggregate flux of distinct through-going line elements.  Sequential recrossing in time remains possible (module177), but simultaneous active-flux divergence in the near-linear catalyst branch cannot be explained by one line piercing the same positive cap repeatedly. '
  'An escape by simultaneous folding must substantially depart from the monotone linear-catalyst geometry and introduce a new non-affine/folding structure; that structure becomes the next object to audit.'),
 'rows':rows
},indent=2,allow_nan=False))
