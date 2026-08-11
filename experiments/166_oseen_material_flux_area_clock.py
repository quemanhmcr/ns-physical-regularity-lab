import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); one=arb(1); rows=[]

# Exact Lamb-Oseen material-circle calibration. Normalize t0=nu=Gamma_infty=1;
# all reported ratios are dimensionless and independent of these choices.
for qs in ('1e-6','1e-3','0.01','0.03','0.1','0.25','0.5','1','2','4','8'):
    q=arb(qs)
    G0=one-(-q).exp()
    # omega(R,t0)=Gamma_inf/(4 pi nu t0) exp(-q)
    # A_Gamma=Gamma_R/omega=4 pi nu t0 (exp(q)-1)
    tau_over_t0=4*pi*(q.exp()-one)
    t1_over_t0=one+tau_over_t0
    q1=q/t1_over_t0
    G1=one-(-q1).exp()
    fraction=G1/G0
    loss=one-fraction
    # Exact material-loop half-life from G(t_half)=G0/2.
    qhalf=-((one+(-q).exp())/2).log()
    thalf_over_t0=q/qhalf
    half_delay_over_t0=thalf_over_t0-one
    half_clocks=half_delay_over_t0/tau_over_t0
    # Exact Kelvin derivative at t0, normalized by G0/tau_A.
    # Gdot=-q exp(-q) at t0 under Gamma_inf=t0=1.
    Gdot0=-q*(-q).exp()
    kelvin_rate_clock=abs(Gdot0)*tau_over_t0/G0
    if not (G0.lower()>0 and tau_over_t0.lower()>0 and fraction.lower()>0 and fraction.upper()<1):
        raise AssertionError(('Oseen positivity',qs,G0,tau_over_t0,fraction))
    rows.append({
        'q0_R2_over_4nut0':qs,
        'initial_material_loop_circulation_fraction_of_total':str(G0),
        'flux_area_clock_tauA_over_t0':str(tau_over_t0),
        'q_after_one_flux_area_clock':str(q1),
        'material_loop_circulation_fraction_remaining_after_one_flux_area_clock':str(fraction),
        'fraction_lost_after_one_flux_area_clock':str(loss),
        'half_life_in_flux_area_clocks':str(half_clocks),
        'initial_Kelvin_change_over_one_flux_area_clock_normalized_by_initial_circulation':str(kelvin_rate_clock),
    })

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
   'Lamb-Oseen is used as an exact Navier-Stokes calibration of local material-circulation persistence.  A circle of radius R is a material loop.  At time t0 let q=R^2/(4 nu t0), Gamma_R=Gamma_infty(1-exp(-q)), and define the physical local flux-area scale A_Gamma=Gamma_R/omega(R,t0)=4 pi nu t0(exp(q)-1). '
   'After one flux-area clock tau_A=A_Gamma/nu, the same material loop has q1=q/[1+4 pi(exp(q)-1)] and circulation fraction [1-exp(-q1)]/[1-exp(-q)].  The half-life is also closed form via q_half=-log[(1+exp(-q))/2]. '
   'For concentrated/core loops q=O(0.1) or larger, the artifact reports order-one local circulation loss on an order-one or smaller number of flux-area clocks.  Extremely diffuse tail loops q->0 are an explicit escape: their relative circulation changes only slowly on this clock. '
   'Thus the relay premise of module165 is supported as a canonical exact NS core calibration but is not universalized to all material tubes.  A general proof must either show that productive maintenance lineages remain in a core-like transverse regime, or separately control the diffuse-tail geometry.'),
 'rows':rows
},indent=2,allow_nan=False))
