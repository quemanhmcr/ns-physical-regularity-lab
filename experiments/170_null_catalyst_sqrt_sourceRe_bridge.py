import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rt3=arb(3).sqrt(); rows=[]
Ccap=arb(4)*pi/(3*rt3); Cbridge=Ccap*arb(7).sqrt()

# B=diag(2,-1,-1), module168 gives C2[self-Euler source]=(A^2 r^2/7) B.
# A productive core with strain amplitude s and radius epsilon has source Re Rs=s epsilon^2/nu.
# To balance a viscous maintenance demand s nu/epsilon^2, set A^2 epsilon^2/7=s nu/epsilon^2.
for nus in ('1e-18','1','1e18'):
    nu=arb(nus)
    for es in ('1e-3','1e-12','1e-40'):
        e=arb(es)
        for Rs in ('1e-12','1','1e6','1e30','1e100'):
            R=arb(Rs)
            s=R*nu/(e*e)
            A2=7*s*nu/(e**4); A=A2.sqrt()
            J=A2*e*e/7
            demand=s*nu/(e*e)
            if not (J/demand).contains(1): raise AssertionError(('maintenance source',nus,es,Rs,J,demand))
            Gamma=Ccap*A*e**3
            ReG=Gamma/nu
            pred=Cbridge*R.sqrt()
            if not (ReG/pred).contains(1): raise AssertionError(('sqrt bridge',nus,es,Rs,ReG,pred))
            Ecore=(arb(4)*pi/105)*A2*e**7
            Epred=(arb(4)*pi/15)*R*nu*nu*e
            if not (Ecore/Epred).contains(1): raise AssertionError(('core energy',nus,es,Rs,Ecore,Epred))
            rows.append({
                'nu':nus,'epsilon':es,'productive_source_Re_s_epsilon2_over_nu':Rs,
                'productive_strain_s':str(s),'required_null_gradient_amplitude_A':str(A),
                'generated_maintenance_source_rate':str(J),
                'positive_cap_catalyst_circulation_at_core':str(Gamma),
                'catalyst_cap_circulation_Re_Gamma_over_nu':str(ReG),
                'closed_sqrt_sourceRe_prediction':str(pred),
                'GammaRe_over_sqrt_sourceRe':str(ReG/R.sqrt()),
                'bridge_constant_4pi_sqrt7_over_3sqrt3':str(Cbridge),
                'catalyst_kinetic_energy_inside_core':str(Ecore),
                'closed_core_energy_4pi_over15_R_nu2_epsilon':str(Epred),
            })

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
   'Restore the actual productive strain amplitude s to the validated linear transaction-null catalyst.  A core of radius epsilon has source Reynolds R_s=s epsilon^2/nu.  Its viscous maintenance demand is proportional to s nu/epsilon^2.  The axisymmetric null catalyst B=diag(2,-1,-1) supplies C2 source A^2 epsilon^2/7, so exact balance requires A^2=7 s nu/epsilon^4. '
   'The catalyst own positive-cap Stokes circulation at the same radius is Gamma_cat=(4 pi/(3 sqrt3)) A epsilon^3.  Eliminating A gives the exact scale-free bridge Gamma_cat/nu=[4 pi sqrt7/(3 sqrt3)] sqrt(R_s). '
   'Therefore the unit-core observation Gamma_cat/nu->0 is not an escape from the singular high-Re branch: whenever the maintained productive source Reynolds diverges, the circulation Reynolds of the cheapest linear-null catalyst also diverges, with square-root law. '
   'The catalyst kinetic energy inside the core is only (4 pi/15) R_s nu^2 epsilon, so this bridge is not an energy contradiction.  It is a direct conversion from productive source Reynolds to material-circulation-scale demand in a concrete transaction-null maintenance mechanism.'),
 'rows':rows
},indent=2,allow_nan=False))
