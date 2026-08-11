import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rt3=arb(3).sqrt(); rows=[]

# Axisymmetric symmetric-traceless null gradient B=diag(2,-1,-1).
# |B|^2=6 and (B^2)_TF=B, so module168 gives Sdot_v(r)=(A^2 r^2/7) B.
# Exact tangent Hodge velocity u=-(A/3)x cross Bx has
# E(B_L)=(2pi/315) A^2 |B|^2 L^7=(4pi/105)A^2 L^7.
# Positive spherical cap n_x>1/sqrt(3) carries actual vorticity flux
# Gamma_cap=(4pi/(3sqrt3)) A R^3.
for ks in ('1','10','100'):
    k=arb(ks)
    for nus in ('1e-12','1','1e12'):
        nu=arb(nus)
        for E0s in ('1e-12','1','1e12'):
            E0=arb(E0s)
            for es in ('1e-3','1e-9','1e-30','1e-90'):
                e=arb(es)
                A2=7*k*nu/(e**4); A=A2.sqrt()
                jcore=A2*e*e/7
                target=k*nu/(e*e)
                if not (jcore/target).contains(1): raise AssertionError(('core source rate',ks,nus,E0s,es,jcore,target))
                L=(arb(15)*E0*e**4/(arb(4)*pi*k*nu))**(arb(1)/7)
                E=(arb(4)*pi/105)*A2*L**7
                if not (E/E0).contains(1): raise AssertionError(('energy horizon',ks,nus,E0s,es,E,E0))
                capfac=arb(4)*pi/(3*rt3)
                Gcore=capfac*A*e**3
                Ghalo=capfac*A*L**3
                Recore=Gcore/nu; Rehalo=Ghalo/nu
                lam=L/e
                rows.append({
                    'maintenance_source_rate_prefactor_k':ks,'nu':nus,'energy_budget_E0':E0s,'epsilon':es,
                    'null_gradient_amplitude_A':str(A),
                    'generated_productive_strain_rate_coefficient_at_core':str(jcore),
                    'energy_horizon_L':str(L),
                    'halo_broadness_lambda_L_over_epsilon':str(lam),
                    'kinetic_energy':str(E),
                    'core_positive_cap_circulation':str(Gcore),
                    'halo_positive_cap_circulation':str(Ghalo),
                    'core_cap_circulation_Re':str(Recore),
                    'halo_cap_circulation_Re':str(Rehalo),
                    'lambda_times_epsilon_power_3over7':str(lam*(e**(arb(3)/7))),
                    'core_Re_divided_by_epsilon':str(Recore/e),
                    'halo_Re_times_epsilon_power_2over7':str(Rehalo*(e**(arb(2)/7))),
                    'halo_over_core_circulation':str(Ghalo/Gcore),
                    'halo_over_core_closed_lambda_cubed':str(lam**3),
                })
                if not ((Ghalo/Gcore)/(lam**3)).contains(1): raise AssertionError(('cap flux scaling',ks,nus,E0s,es))

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
   'Scale the exact linear transaction-null catalyst of module168 using B=diag(2,-1,-1), for which (B^2)_TF=B.  To generate a productive strain-rate coefficient k nu/epsilon^2 at the maintained core requires A^2=7 k nu/epsilon^4. '
   'The exact tangent Hodge velocity has kinetic energy E=(4 pi/105)A^2 L^7.  Saturating a finite energy budget gives L~epsilon^(4/7), so the null catalyst may be much broader than the core: L/epsilon~epsilon^-3/7. '
   'This geometry has an actual Stokes circulation observable.  On a sphere of radius R the positive cap n_x>1/sqrt(3) carries Gamma_cap=(4 pi/(3 sqrt(3))) A R^3.  At the tiny core Gamma_core/nu~epsilon ->0, but at the finite-energy halo horizon Gamma_halo/nu~epsilon^-2/7 ->infinity. '
   'Thus the linear null catalyst genuinely evades the need for pre-existing productive transaction at the core, but it does not remain a low-circulation ancestry mechanism when broadened enough to exploit the Oseen reservoir escape.  It relocates circulation ancestry from the tiny core into a shrinking but increasingly high-Re halo. '
   'This is a finite-energy calibration, not yet a theorem that the cap flux is one closed frozen lineage.  The next attack must compare this required halo flux with the exact closed-lineage specific-volume bound; if it exceeds that capacity, the catalyst is forced into through-going recruitment or viscous renewal.'),
 'rows':rows
},indent=2,allow_nan=False))
