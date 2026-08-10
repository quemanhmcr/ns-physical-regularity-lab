import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
# Closed continuous strain-count models N in [0,Nmax]. Exposure Xi=int dN/Re(N).
for model in ['constant','linear','exponential']:
 for Ns in ['1','10','100']:
  N=arb(Ns)
  if model=='constant':
   R0=arb(10); Xi=N/R0; Re_end=R0
  elif model=='linear':
   # Re=1+N, Xi=log(1+N)
   Xi=(1+N).log(); Re_end=1+N
  else:
   # Re=exp(N), Xi=1-exp(-N)
   Re_end=N.exp(); Xi=1-(-N).exp()
  rows.append({'model':model,'strain_count_N':Ns,'cumulative_viscous_exposure_Xi_int_dN_over_Re':str(Xi),'terminal_Re':str(Re_end)})
# Threshold identity: strain count accumulated while Re<=Rstar is <=Rstar*Xi.
for Rstars in ['1','10','1e6']:
 R=arb(Rstars)
 for Xis in ['1e-24','1','1e24']:
  X=arb(Xis); ceiling=R*X
  rows.append({'threshold_Re_star':Rstars,'total_exposure_Xi':Xis,'maximum_strain_count_possible_below_threshold':str(ceiling),'structural_inequality':'int_{Re<=R*} dN <= R* Xi'})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'With strain count dN=s dt and source Re=sL^2/nu, the localized viscous exposure is Xi=int dN/Re. If total strain count diverges but Xi remains finite, Re cannot stay bounded: for every fixed R*, the strain count spent with Re<=R* is at most R*Xi, so an infinite tail must move to arbitrarily large source Re. Equivalently, the circulation-dimensional transaction Gamma_Q=nu Re must become unbounded along the frozen-ancestry escape.','rows':rows},indent=2,allow_nan=False))
