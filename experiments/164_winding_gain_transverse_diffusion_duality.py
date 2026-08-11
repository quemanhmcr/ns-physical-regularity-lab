import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rt2=arb(2).sqrt(); rows=[]; target_rows=[]

# Exact spherical-shell winding field from modules 104-105.
# omega_r=A sin^2(theta)cos(theta)/r^2, omega_phi=B sin(theta)cos(theta)/r.
# At theta=pi/4: omega_r=A/(2 sqrt(2) r^2), omega_phi=B/(2r).
# Positive north-hemisphere throughgoing circulation is Gamma=pi A/2.
# N turns across shell DeltaR gives B/A=pi N/DeltaR.
# Hodge transaction circulation is Gamma_Q=2Br/5.
for R1s,R2s in [('1','2'),('1e-6','2e-6'),('1e6','2e6')]:
    R1=arb(R1s); R2=arb(R2s); dR=R2-R1; r=(R1*R2).sqrt()
    for Gs in ('1e-18','1','1e18'):
        Gamma=arb(Gs); A=2*Gamma/pi
        for nus in ('1e-12','1','1e12'):
            nu=arb(nus)
            for Ns in ('1','10','1e4','1e12'):
                N=arb(Ns); B=A*pi*N/dR
                GammaQ=2*B*r/5
                gain=GammaQ/Gamma
                omr=A/(2*rt2*r*r); omp=B/(2*r); om=(omr*omr+omp*omp).sqrt()
                # Characteristic transverse flux-area scale at the representative point.
                Aflux=abs(Gamma)/om
                tau_nu=Aflux/nu
                tau_trans=r*r/abs(GammaQ)
                persistence=tau_nu/tau_trans
                correction=(1+dR*dR/(2*pi*pi*N*N*r*r)).sqrt()
                closed=(arb(4)/5)*(abs(Gamma)/nu)/correction
                if not (persistence/closed).contains(1):
                    raise AssertionError(('winding persistence identity',R1s,R2s,Gs,nus,Ns,persistence,closed))
                limit=(arb(4)/5)*(abs(Gamma)/nu)
                null_ratio=(arb(32)*pi/105*A*A*dR/(R1*R2))/(arb(8)*pi/15*B*B*dR)
                rows.append({
                    'R1':R1s,'R2':R2s,'Gamma':Gs,'nu':nus,'turns_N':Ns,
                    'transaction_gain_GammaQ_over_Gamma':str(gain),
                    'lineage_circulation_Re_Gamma_over_nu':str(abs(Gamma)/nu),
                    'characteristic_flux_area_Gamma_over_abs_omega':str(Aflux),
                    'viscous_clock_area_over_nu':str(tau_nu),
                    'transaction_clock_r2_over_GammaQ':str(tau_trans),
                    'viscous_over_transaction_clock':str(persistence),
                    'closed_clock_ratio':str(closed),
                    'high_winding_limit_4Gamma_over_5nu':str(limit),
                    'clock_ratio_over_high_winding_limit':str(persistence/limit),
                    'transaction_null_over_sharp_enstrophy':str(null_ratio),
                })

# Adversarial form: hold target transaction Reynolds fixed while lineage circulation Reynolds ->0.
# Solve N exactly from Gamma_Q/nu = RQ using DeltaR/r geometry.
R1=arb(1);R2=arb(2);dR=R2-R1;r=(R1*R2).sqrt();nu=arb(1)
for RQs in ('1','10','30'):
    RQ=arb(RQs)
    for RGs in ('1','1e-2','1e-6','1e-12','1e-30'):
        RG=arb(RGs);Gamma=RG*nu;A=2*Gamma/pi
        N=(arb(5)*dR/(4*r))*(RQ/RG)
        B=A*pi*N/dR;GammaQ=2*B*r/5
        if not (GammaQ/(nu*RQ)).contains(1): raise AssertionError(('target transaction',RQs,RGs,GammaQ))
        omr=A/(2*rt2*r*r);omp=B/(2*r);om=(omr*omr+omp*omp).sqrt()
        Aflux=Gamma/om;tau_nu=Aflux/nu;tau_trans=r*r/GammaQ;persistence=tau_nu/tau_trans
        scaled=persistence/RG
        high=arb(4)/5
        # scaled tends 4/5 as RG/RQ ->0; report exact finite correction.
        target_rows.append({
            'target_transaction_Re_GammaQ_over_nu':RQs,
            'lineage_circulation_Re_Gamma_over_nu':RGs,
            'required_turns_N':str(N),
            'viscous_over_transaction_clock':str(persistence),
            'clock_ratio_divided_by_lineage_Re':str(scaled),
            'high_winding_limit_of_scaled_ratio_4_over_5':str(high),
            'scaled_ratio_over_4_over_5':str(scaled/high),
        })

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'target_cases':len(target_rows),
 'interpretation':(
   'The exact slow-spiral shell is used to attack the inference that large winding gain can let a vanishing-circulation material lineage continuously finance an order-one maintenance transaction. '
   'At theta=pi/4 the field has |omega|^2=B^2/(4r^2)+A^2/(8r^4).  The throughgoing circulation Gamma defines the characteristic transverse flux-area scale A_Gamma=Gamma/|omega|, whose viscous clock is A_Gamma/nu.  The Hodge transaction clock is r^2/Gamma_Q. '
   'Their exact ratio is (4/5)(Gamma/nu)[1+DeltaR^2/(2 pi^2 N^2 r^2)]^-1/2.  Thus as winding N grows, Gamma_Q/Gamma grows linearly and the transaction-null fraction falls quadratically, but the viscous-to-transaction clock ratio approaches only (4/5)Gamma/nu and does not inherit the winding gain. '
   'Equivalently, at fixed target transaction Reynolds Gamma_Q/nu, choosing lineage circulation Reynolds Gamma/nu ->0 forces N->infinity while the characteristic persistence ratio tends to zero linearly with Gamma/nu. '
   'This is an exact calibration identity for the validated spherical-shell winding geometry, not yet a universal theorem for arbitrary flux-tube shapes.  It kills the simplest vanishing-circulation winding escape for continuous maintenance over one transaction clock.  The surviving escapes require continual lineage relay/viscous ancestry renewal, geometrically nonuniform tubes that evade this calibration, or a different nonlinear maintenance mechanism.'),
 'rows':rows,'fixed_transaction_rows':target_rows
},indent=2,allow_nan=False))
