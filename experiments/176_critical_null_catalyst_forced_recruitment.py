import json, os
from fractions import Fraction as F
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rt3=arb(3).sqrt()
# Fixed p=2 turnover used in modules172-173.
IE=F(1,7)-F(4,9)+F(6,11)-F(4,13)+F(1,15)
def A(q):return arb(q.numerator)/q.denominator
CE=(arb(4)*pi/15)*A(IE)
speak=(arb(3)/7).sqrt(); fpeak=(speak**3)*(1-speak*speak)**2
Ccap=arb(4)*pi/(3*rt3); kappa=arb(77)/2
rows=[]
for alphas in ('0.4','0.45','0.49'):
 alpha=arb(alphas); beta=(1+4*alpha)/7
 for nus in ('1e-12','1','1e12'):
  nu=arb(nus)
  for E0s in ('1e-6','1','1e6'):
   E0=arb(E0s)
   # tau=1 reference for the same power-law branch.
   tau0=arb(1);eps0=tau0**alpha;A20=7*nu/(tau0*eps0**4);L0=(E0/(CE*A20))**(arb(1)/7);G0=Ccap*A20.sqrt()*L0**3*fpeak
   for taus in ('1e-6','1e-30','1e-100'):
    tau=arb(taus);eps=tau**alpha;strain=1/tau;A2=7*strain*nu/(eps**4);Agr=A2.sqrt();L=(E0/(CE*A2))**(arb(1)/7)
    Gamma=Ccap*Agr*L**3*fpeak
    direct=(beta/2)*Gamma/tau
    shrink=-3*beta*Gamma/tau
    kelvin=-kappa*nu*Gamma/(L*L)
    selfe=arb(0)
    recruit=direct-selfe-shrink-kelvin
    closednorm=7*beta/2+kappa*nu*tau/(L*L)
    if not ((recruit/(Gamma/tau))/closednorm).contains(1): raise AssertionError(('recruit identity',alphas,nus,E0s,taus))
    lower_cumulative=7*(Gamma-G0)
    if tau<1 and not (lower_cumulative.lower()>0): raise AssertionError(('positive cumulative lower',alphas,nus,E0s,taus,lower_cumulative))
    rows.append({
      'alpha':alphas,'beta':str(beta),'nu':nus,'E0':E0s,'tau':taus,'epsilon':str(eps),'L':str(L),
      'peak_cap_circulation':str(Gamma),'direct_required_cap_flux_growth_rate':str(direct),
      'self_Euler_flux_rate':'0','geometric_shrinkage_recruitment_rate':str(shrink),'turnover_Kelvin_flux_rate':str(kelvin),
      'required_external_source_relative_recruitment_rate':str(recruit),
      'external_rate_divided_by_Gamma_over_tau':str(recruit/(Gamma/tau)),
      'asymptotic_leading_coefficient_7beta_over2':str(7*beta/2),
      'Kelvin_correction_coefficient_kappa_nu_tau_over_L2':str(kappa*nu*tau/(L*L)),
      'Kelvin_magnitude_over_direct_growth':str(abs(kelvin)/direct),
      'external_rate_over_direct_growth':str(recruit/direct),
      'cumulative_external_recruitment_lower_bound_from_tau1':str(lower_cumulative),
      'lower_bound_over_current_Gamma':str(lower_cumulative/Gamma),
    })
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'Combine the critical null-catalyst scaling with the exact source-relative moving-cap and Kelvin ledgers.  The peak cap circulation obeys Gamma~tau^(-beta/2), beta=(1+4alpha)/7, so its required direct growth is +(beta/2)Gamma/tau. '
  'The catalyst self-Euler conversion changes cap flux by zero.  The shrinking active cap contributes -3 beta Gamma/tau, and the p=2 turnover collar contributes the negative Kelvin rate -(77/2)nu Gamma/L^2.  Therefore the missing positive source-relative material recruitment is exactly [7 beta/2 +(77/2)nu tau/L^2] Gamma/tau. '
  'Because beta<1/2, nu tau/L^2->0 in the critical branch, so recruitment asymptotically approaches (7 beta/2)Gamma/tau and is seven times the net cap-flux growth rate.  Integrating only this leading positive requirement from tau=1 gives the rigorous branchwise lower bound 7[Gamma(tau)-Gamma(1)], which diverges with Gamma. '
  'This is the first infinite required ancestry-throughput statement produced by the physical moving-source ledger rather than by recurrence counting.  It is not yet a contradiction: the same material circulation may in principle cross the moving source repeatedly.  The next attack must therefore test recross/reuse rather than search for another scalar budget.'),
 'rows':rows
},indent=2,allow_nan=False))
