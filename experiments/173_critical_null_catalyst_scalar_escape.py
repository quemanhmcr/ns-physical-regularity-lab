import json, os
from fractions import Fraction as F
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rt3=arb(3).sqrt()
# Fixed localization chi=(1-s^2)^2 from module172.
IE=F(1,7)-F(4,9)+F(6,11)-F(4,13)+F(1,15)
ID=16*(F(1,9)-F(2,11)+F(1,13))
def A(q):return arb(q.numerator)/q.denominator
CE=(arb(4)*pi/15)*A(IE)
CZ=(arb(8)*pi/15)*A(ID)
speak=(arb(3)/arb(7)).sqrt(); fpeak=(speak**3)*(1-speak*speak)**2
Ccap=arb(4)*pi/(3*rt3)
rows=[]
for alphas in ('0.4','0.45','0.49'):
    alpha=arb(alphas); beta=(1+4*alpha)/7
    if not (beta<arb('0.5')): raise AssertionError(('beta',alphas,beta))
    tail_exp=1-2*beta
    for nus in ('1e-12','1','1e12'):
        nu=arb(nus)
        for E0s in ('1e-6','1','1e6'):
            E0=arb(E0s)
            for taus in ('1e-6','1e-30','1e-100'):
                tau=arb(taus); eps=tau**alpha; strain=1/tau
                A2=7*strain*nu/(eps**4); Agr=A2.sqrt()
                L=(E0/(CE*A2))**(arb(1)/7)
                E=CE*A2*L**7; Z=CZ*A2*L**5
                if not (E/E0).contains(1): raise AssertionError(('energy',alphas,nus,E0s,taus,E,E0))
                lifetime=E/(nu*Z)
                # Z(t)=const*tau^(-2 beta), so tail integral from 0 to tau is nu Z(tau) tau/(1-2beta).
                diss_tail=nu*Z*tau/tail_exp
                Gamma_peak=Ccap*Agr*L**3*fpeak; ReGamma=Gamma_peak/nu
                Rs=strain*eps*eps/nu
                rows.append({
                    'alpha_core_radius_exponent':alphas,'beta_halo_radius_exponent':str(beta),'nu':nus,'energy_budget_E0':E0s,'tau':taus,
                    'epsilon':str(eps),'strain_s':str(strain),'source_Re_s_epsilon2_over_nu':str(Rs),
                    'null_gradient_A':str(Agr),'energy_horizon_L':str(L),'halo_over_core_L_over_epsilon':str(L/eps),
                    'catalyst_energy':str(E),'turnover_enstrophy':str(Z),'remaining_viscous_dissipation_tail':str(diss_tail),
                    'tail_exponent_1_minus_2beta':str(tail_exp),'halo_energy_dissipation_clock':str(lifetime),
                    'halo_clock_over_remaining_time':str(lifetime/tau),'peak_positive_cap_circulation_Re':str(ReGamma),
                    'ReGamma_times_tau_power_beta_over2':str(ReGamma*(tau**(beta/2))),
                    'L_divided_by_tau_power_beta':str(L/(tau**beta)),'Z_times_tau_power_2beta':str(Z*(tau**(2*beta))),
                    'dissipation_tail_divided_by_tau_power_1minus2beta':str(diss_tail/(tau**tail_exp)),
                })
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'fixed_taper_energy_constant_CE':str(CE),'fixed_taper_enstrophy_constant_CZ':str(CZ),
 'interpretation':(
  'Insert the singular scaling s=1/tau and epsilon=tau^alpha, 2/5<=alpha<1/2, into the validated linear null catalyst and localize it with the fixed physical velocity taper chi=(1-r^2/L^2)^2.  The required gradient obeys A^2=7 nu tau^(-1-4alpha). '
  'Holding the catalyst kinetic occupancy at finite E0 gives L~tau^beta with beta=(1+4alpha)/7<1/2.  The exact turnover identity then gives Z~tau^(-2beta), so the remaining viscous dissipation tail nu int_0^tau Z dt is finite and vanishes like tau^(1-2beta). '
  'Meanwhile the catalyst energy-dissipation clock E/(nu Z)~L^2/nu has ratio to the remaining singular time proportional to tau^(2beta-1), which diverges.  The localized halo therefore becomes effectively frozen over the time remaining to the putative singularity even though its closure collar has nonzero viscosity. '
  'The positive-cap circulation Reynolds of the halo diverges as tau^(-beta/2).  Thus turnover viscosity does not provide the missing scalar contradiction in the admissible source corridor; it leaves an increasingly high-circulation through-going ancestry reservoir whose assembly/recruitment is the unresolved burden.  This is an adversarial scaling calibration, not an exact blow-up solution.'),
 'rows':rows
},indent=2,allow_nan=False))
