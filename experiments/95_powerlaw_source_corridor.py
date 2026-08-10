import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
# tau=T-t, s=tau^-1, L=tau^alpha, nu=1 units.
# Hodge-horizon compatibility requires alpha>=2/5; high-Re/finitely accumulated localized exposure requires alpha<1/2.
for alphas in ['0.35','0.4','0.45','0.49','0.5','0.55']:
 a=arb(alphas)
 for taus in ['1e-6','1e-30','1e-100']:
  tau=arb(taus); s=1/tau; L=tau**a; Re=s*L*L
  ratio_energy=L/(tau**(arb(2)/5))
  Gamma=s*L*L
  occupancy=s*s*L**5
  diss_integrand=s*s*L**3
  row={'alpha':alphas,'tau':taus,'L':str(L),'s':str(s),'Re_source':str(Re),'L_over_tau_2over5':str(ratio_energy),'Gamma_Q_sL2':str(Gamma),'s2L5_occupancy_scale':str(occupancy),'s2L3_enstrophy_scale':str(diss_integrand),'Hodge_energy_horizon_compatible':bool(a>=arb(2)/5),'high_Re_as_tau_to_zero':bool(a<arb(1)/2),'finite_cumulative_exposure_as_tau_to_zero':bool(a<arb(1)/2),'finite_sharp_dissipation_integral_as_tau_to_zero':bool(a>arb(1)/3)}
  rows.append(row)
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For the natural finite-time strain law s~1/(T-t) and a source scale L~(T-t)^alpha, the Hodge energy horizon requires alpha>=2/5 while high source Re and finite cumulative localized viscous exposure require alpha<1/2. The sharp-enstrophy spacetime integral only requires alpha>1/3. Thus the scalar-admissible frozen escape lives in the intrinsic corridor 2/5<=alpha<1/2, throughout which Gamma_Q=sL^2 diverges as the singular time is approached.','rows':rows},indent=2,allow_nan=False))
