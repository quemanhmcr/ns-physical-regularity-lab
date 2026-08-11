import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rt2=arb(2).sqrt(); rows=[]

# Same exact shell calibration as module164.
# Conditional relay throughput: if one lineage is trusted for at most tau_nu,
# the number of lineage-lifetimes required to cover one transaction clock is tau_trans/tau_nu.
# Multiplying by the circulation per lineage gives A_relay = Gamma tau_trans/tau_nu.
for R1s,R2s in [('1','2'),('1e-6','2e-6'),('1e6','2e6')]:
    R1=arb(R1s);R2=arb(R2s);dR=R2-R1;r=(R1*R2).sqrt()
    for nus in ('1e-18','1','1e18'):
        nu=arb(nus)
        for RQs in ('1','10','30'):
            RQ=arb(RQs)
            for RGs in ('1','1e-2','1e-6','1e-12','1e-30'):
                RG=arb(RGs);Gamma=RG*nu
                N=(arb(5)*dR/(4*r))*(RQ/RG)
                A=2*Gamma/pi; B=A*pi*N/dR; GammaQ=2*B*r/5
                if not (GammaQ/(nu*RQ)).contains(1): raise AssertionError(('target transaction',R1s,R2s,nus,RQs,RGs))
                omr=A/(2*rt2*r*r);omp=B/(2*r);om=(omr*omr+omp*omp).sqrt()
                Aflux=abs(Gamma)/om;tau_nu=Aflux/nu;tau_trans=r*r/abs(GammaQ)
                persist=tau_nu/tau_trans
                lifetimes=tau_trans/tau_nu
                relay=abs(Gamma)*lifetimes
                correction=(1+dR*dR/(2*pi*pi*N*N*r*r)).sqrt()
                closed=(arb(5)/4)*nu*correction
                if not (relay/closed).contains(1): raise AssertionError(('relay identity',R1s,R2s,nus,RQs,RGs,relay,closed))
                rows.append({
                    'R1':R1s,'R2':R2s,'nu':nus,
                    'target_transaction_Re_GammaQ_over_nu':RQs,
                    'lineage_circulation_Re_Gamma_over_nu':RGs,
                    'required_winding_turns_N':str(N),
                    'viscous_over_transaction_clock':str(persist),
                    'lineage_lifetimes_per_transaction_clock':str(lifetimes),
                    'conditional_relay_circulation_throughput':str(relay),
                    'relay_throughput_over_nu':str(relay/nu),
                    'closed_relay_throughput_over_nu':str(closed/nu),
                    'high_winding_limit_relay_throughput_over_nu':'1.25',
                })

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
   'Use the exact winding/diffusion identity from module164 and ask a conditional spacetime question: if one throughgoing lineage can be trusted for at most one characteristic transverse viscous clock tau_nu, how much circulation ancestry must be presented by successive lineage-lifetimes to cover one transaction clock tau_trans? '
   'The exact shell answer is A_relay=Gamma tau_trans/tau_nu=(5/4) nu sqrt[1+DeltaR^2/(2 pi^2 N^2 r^2)].  Thus at fixed target transaction Reynolds, taking Gamma/nu->0 drives N and the required number of lineage-lifetimes to infinity, but their aggregate circulation throughput per transaction clock tends to (5/4)nu rather than zero. '
   'This quantity is not promoted as an irreversible energy or circulation cost.  It is a conditional ancestry-throughput demand.  In the exact moving-source material-spacetime ledger, such replacement can only be served by relative material crossing of the source boundary or by Kelvin viscous current. '
   'The premise that one transverse diffusion clock really destroys an order-one fraction of useful local material circulation must be attacked separately in an exact Navier-Stokes solution; the next module uses Lamb-Oseen for that calibration.'),
 'rows':rows
},indent=2,allow_nan=False))
