import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rows=[]

# Capacity stationary productive ray S_*=sigma S_unit with |S_unit|^2=21/2.
Sunit2=arb(21)/2
# Exact tangent carrier: S_v(r)=(r/L)^2 S_outer,
# E=(8pi/405)|S_outer|^2 L^5, Z=(112pi/45)|S_outer|^2 L^3.
for cs in ('1','10','42'):
    c=arb(cs)
    for nus in ('1e-12','1','1e12'):
        nu=arb(nus)
        for E0s in ('1e-12','1','1e12'):
            E0=arb(E0s)
            for es in ('1e-3','1e-9','1e-30','1e-90'):
                e=arb(es)
                sigma_core=c*nu/(e*e)
                # Energy-horizon broadness from E=(28pi/135)c^2 nu^2 e lambda^9.
                lam_max=(arb(135)*E0/(arb(28)*pi*c*c*nu*nu*e))**(arb(1)/9)
                Lmax=e*lam_max
                sigma_outer=sigma_core*lam_max*lam_max
                Re_outer=sigma_outer*Lmax*Lmax/nu
                E=(arb(8)*pi/405)*(Sunit2*sigma_outer*sigma_outer)*Lmax**5
                Z=(arb(112)*pi/45)*(Sunit2*sigma_outer*sigma_outer)*Lmax**3
                if not (E/E0).contains(1): raise AssertionError(('energy horizon',cs,nus,E0s,es,E,E0))
                closed_Re=c*lam_max**4
                if not (Re_outer/closed_Re).contains(1): raise AssertionError(('outer Re',cs,nus,E0s,es,Re_outer,closed_Re))
                q_like=(e/Lmax)**2
                rows.append({
                    'maintenance_gate_c_sigma_e2_over_nu':cs,'nu':nus,'energy_budget_E0':E0s,'epsilon':es,
                    'energy_horizon_lambda_L_over_epsilon':str(lam_max),
                    'energy_horizon_L':str(Lmax),
                    'local_maintenance_sigma':str(sigma_core),
                    'outer_carrier_sigma_needed':str(sigma_outer),
                    'outer_productive_transaction_Re_sigma_outer_L2_over_nu':str(Re_outer),
                    'closed_outer_Re_c_lambda4':str(closed_Re),
                    'carrier_kinetic_energy':str(E),
                    'carrier_enstrophy':str(Z),
                    'deep_core_parameter_q_like_epsilon2_over_L2':str(q_like),
                    'lambda_times_epsilon_power_1over9':str(lam_max*(e**(arb(1)/9))),
                    'outer_Re_times_epsilon_power_4over9':str(Re_outer*(e**(arb(4)/9))),
                    'q_like_divided_by_epsilon_power_2over9':str(q_like/(e**(arb(2)/9))),
                })

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
   'Attack the broad-circulation-reservoir escape with the cheapest smooth productive l=2 tangent carrier.  Let epsilon be the maintained core radius and L=lambda epsilon the smooth productive halo radius.  Because the exact carrier has S_v(epsilon)=(epsilon/L)^2 S_outer, producing the maintenance strain sigma_core=c nu/epsilon^2 requires sigma_outer=lambda^2 sigma_core. '
   'For the stationary productive ray |S_unit|^2=21/2, the exact tangent-carrier energy becomes E=(28 pi/135)c^2 nu^2 epsilon lambda^9, while the halo transaction Reynolds is R_halo=sigma_outer L^2/nu=c lambda^4. '
   'Thus a core-like reservoir lambda=O(1) stays order-one Re and is exposed to the material-circulation renewal calibrated by Oseen/Burgers.  A broad reservoir lambda->infinity can make the local ancestry parameter q_like=(epsilon/L)^2 small, but necessarily drives its own productive transaction Reynolds to infinity like lambda^4. '
   'Finite kinetic energy does not kill this escape: at the energy horizon lambda~epsilon^-1/9, q_like~epsilon^2/9 and R_halo~epsilon^-4/9.  Instead it routes the broad-reservoir escape back into the previously isolated high-Re shrinking-source ancestry branch. '
   'This is an exact calibration for the minimum smooth tangent carrier, not a universal theorem that every broad vorticity halo has this energy.  Extra transaction-null structure can only add angular enstrophy at fixed Q but may alter kinetic-energy geometry, so the next attack must keep the distinction between productive halo and arbitrary ancestry reservoir.'),
 'rows':rows
},indent=2,allow_nan=False))
