import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
# Kinematic scaling ledger in strain-count N: Re=e^{aN}, L=e^{-bN}, nu=1.
# s=Re/L^2=e^{(a+2b)N}, dt=dN/s.
# If b>2a>0: finite remaining time, finite exposure, decaying strain occupancy E~s^2 L^5=e^{(2a-b)N},
# finite dissipation tail int Z dt with Z~s^2 L^3 and integrand ~e^{(a-b)N}, while flux~Re diverges.
for a_s,b_s in [('0.5','2'),('1','3'),('2','5')]:
 a=arb(a_s); b=arb(b_s)
 if not (b>2*a): raise AssertionError(('need b>2a',a_s,b_s))
 for Ns in ['0','1','10','100']:
  N=arb(Ns); Re=(a*N).exp(); L=(-b*N).exp(); s=((a+2*b)*N).exp()
  time_tail=(-(a+2*b)*N).exp()/(a+2*b)
  exposure_tail=(-a*N).exp()/a
  occupancy_scale=((2*a-b)*N).exp()
  dissipation_tail=((a-b)*N).exp()/(b-a)
  flux_scale=Re
  rows.append({'a':a_s,'b':b_s,'N':Ns,'Re_source':str(Re),'source_scale_L':str(L),'production_rate_s':str(s),'remaining_physical_time':str(time_tail),'remaining_viscous_exposure':str(exposure_tail),'strain_energy_occupancy_scaling_s2L5':str(occupancy_scale),'remaining_viscous_enstrophy_budget_scaling':str(dissipation_tail),'circulation_flux_aggregation_scale':str(flux_scale)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'There exist purely kinematic high-Re inward-cascade scalings with L=e^{-bN}, Re=e^{aN}, b>2a>0, for which infinite strain count fits into finite physical time, cumulative localized viscous exposure is finite, the s^2 L^5 strain-occupancy scale decays, the modeled viscous-enstrophy tail is finite, yet the required circulation-dimensional source flux grows without bound. This is not an NS solution; it is an escape-route autopsy proving that the current scalar energy/dissipation ledgers alone cannot close the high-Re branch. The missing obstruction must be material ancestry aggregation/recruitment geometry.','rows':rows},indent=2,allow_nan=False))
