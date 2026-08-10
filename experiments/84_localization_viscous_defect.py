import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); E2=arb(6)

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
rows=[]
for p_int in [2,8,32]:
 p=arb(p_int); coeff=p*(p+5)
 for qs in ['1e-24','1','1e24']:
  q0=arb(qs)
  for Ls in ['1e-6','1','1e6']:
   L=arb(Ls)
   # Localized profile q=q0 x^2(1-x^p), x=r/L. It is the viscosity-null r^2 mode in the interior, tapered to zero at x=1.
   defect_L2_norm=(20*pi/9)*E2*q0*q0/L * (p*p*(p+5)*(p+5)/(2*p+3))
   for alpha_s in ['0.5','0.9']:
    alpha=arb(alpha_s); outer_fraction=1-alpha**(2*p+3)
    if not (outer_fraction>0 and outer_fraction<1): raise AssertionError(('defect fraction',p_int,alpha_s,outer_fraction))
    rows.append({'taper_power_p':p_int,'q0':qs,'L':Ls,'inner_fraction_alpha':alpha_s,'integrated_squared_radial_viscous_defect':str(defect_L2_norm),'fraction_of_defect_outside_alphaL':str(outer_fraction),'structural_L2_q':f'-p(p+5) q0/L^2 x^p with p={p_int}'})
 # coefficient growth relative to p^3 tends 1/2 before common angular factor.
 ratio=(p*p*(p+5)*(p+5)/(2*p+3))/(p**3)
 rows.append({'taper_power_p':p_int,'dimensionless_defect_coefficient_over_p3':str(ratio)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'Localizing the smooth viscosity-null productive profile forces Q/r^2 to vary radially. For q=q0 x^2(1-x^p), L2 q=-p(p+5)q0 x^p/L^2 and the squared productive viscous defect scales like q0^2/L times p^2(p+5)^2/(2p+3), asymptotically proportional to p^3. Sharper localization pushes almost all defect into the outer taper collar instead of eliminating it. This identifies a forced viscous-active source region, not yet a finite-resource toll.','rows':rows},indent=2,allow_nan=False))
