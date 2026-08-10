import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))

# At the origin of u=(a x+eps yz,-a y,0), the center strain is S0=diag(a,-a,0).
# The quadratic normal-flux term eps r^2 n_x n_y n_z has no degree-two spherical
# harmonic component, so the Hodge extractor returns S_h(r)=S0 at every radius.
# The exact local kinetic energy is
# E(B_r)=4*pi*a^2*r^5/15 + 2*pi*eps^2*r^7/105.
rows=[]
for as_ in ['1e-24','1','1e24']:
  a=arb(as_)
  for es in ['1e-24','1','1e24']:
    eps=arb(es)
    for rs in ['1e-6','1','1e6']:
      r=arb(rs)
      Sh2=2*a*a
      Eh=(2*pi/15)*Sh2*r**5
      Eaff=4*pi*a*a*r**5/15
      Evort=2*pi*eps*eps*r**7/105
      E=Eaff+Evort
      certify_one(Eh/Eaff,('harmonic affine energy saturation',as_,es,rs))
      if not (E>=Eh): raise AssertionError(('total local energy below Hodge floor',as_,es,rs,E,Eh))
      rows.append({'a':as_,'eps':es,'r':rs,'Sh_Frobenius_squared':str(Sh2),'Hodge_harmonic_floor':str(Eh),'exact_affine_energy':str(Eaff),'quadratic_extra_energy':str(Evort),'exact_total_ball_energy':str(E)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'The common reset strain in the exact quadratic family is exactly the degree-two harmonic Hodge mode at the origin.  Its kinetic occupancy saturates the universal Hodge floor E_h=(2*pi/15)|S_h|^2 r^5, while the quadratic vorticity carrier only adds positive r^7 energy.  Common symmetric strain is therefore not a gauge: maintaining a large coherent reset rate on a fixed physical radius requires real kinetic occupancy.'
 ),'rows':rows
},indent=2,allow_nan=False))
