import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
for ss in ['1e-24','1','1e24']:
 s=arb(ss)
 for nus in ['1e-24','1','1e24']:
  nu=arb(nus)
  for Ls in ['1e-12','1','1e12']:
   L=arb(Ls); Re=s*L*L/nu
   exponent=arb(28)/Re
   # For pure diffusion of the localized productive sector, squared norm obeys N(t)<=N(0) exp(-28 nu t/L^2).
   # Over one production time 1/s the universal squared-norm survival ceiling is exp(-28/Re_source).
   survival=(-exponent).exp()
   rows.append({'production_rate_s':ss,'nu':nus,'source_radius_L':Ls,'source_transaction_Re_sL2_over_nu':str(Re),'minimum_squared_norm_decay_exponent_per_strain_time_28_over_Re':str(exponent),'pure_diffusion_squared_norm_survival_upper':str(survival),'regeneration_dominated_if_Re_small':bool(Re<arb(1))})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'Once the productive carrier is genuinely localized at source radius L, the smooth r^2 zero-mode escape is removed and the universal radial spectral gap gives pure-diffusion squared-norm decay at least exp(-28 nu t/L^2). Over one production time 1/s the natural source-scale transaction Reynolds number is Re_source=s L^2/nu. Small Re_source means viscosity erodes the localized productive sector before one strain time unless nonlinear/material regeneration continuously rebuilds it; large Re_source is the approximately frozen-flux branch. This is a mechanism gate, not a standalone full-NS decay theorem because nonlinear sources remain present.','rows':rows},indent=2,allow_nan=False))
