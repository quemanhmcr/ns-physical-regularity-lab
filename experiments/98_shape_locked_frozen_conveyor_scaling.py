import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
# Conditional shape-lock/no-renewal kinematics: d log rho_a=d log rho_b=dN, d log r=-2dN.
for alphas in ['0.4','0.45','0.49']:
 alpha=arb(alphas)
 for taus in ['1e-6','1e-30','1e-100']:
  tau=arb(taus); N=-(alpha/2)*tau.log(); rho=N.exp(); r=(-2*N).exp(); s=alpha/(2*tau); Re=s*r*r # nu=1 units
  if not (r/(tau**alpha)).contains(1): raise AssertionError(('bridge power law',alphas,taus,r))
  if not ((rho*rho*r).contains(1)): raise AssertionError(('oriented cell magnitude preserved',alphas,taus,rho,r))
  rows.append({'alpha_bridge':alphas,'tau':taus,'strain_count_N':str(N),'endpoint_vorticity_amplification_rho_over_rho0':str(rho),'bridge_ratio_r_over_r0':str(r),'rho_squared_times_r':str(rho*rho*r),'shape_locked_gain_rate_s_alpha_over_2tau':str(s),'bridge_source_Re_nu1':str(Re),'endpoint_blowup_exponent_alpha_over_2':str(alpha/2)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'Conditionally on maintaining productive Gram shape and negligible pair-cell renewal, equal endpoint amplification dN forces bridge contraction e^{-2N}. If the bridge/source scale follows r~tau^alpha, then each endpoint magnitude grows like tau^{-alpha/2} and the gain rate is alpha/(2tau). This is a kinematic frozen-conveyor calibration, not an exact finite-energy NS blow-up solution; it exposes how incompressibility can trade two stretching directions against one shrinking source direction without changing the oriented pair-cell determinant.','rows':rows},indent=2,allow_nan=False))
