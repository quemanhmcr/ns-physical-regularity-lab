import json, os, math
import mpmath as mp

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
mp.mp.dps=int(BITS*math.log10(2))+35
pi=mp.pi

d=mp.mpf('1')
I0=mp.mpf('1')            # fixed impulse magnitude of each ring carrier
core_fraction=mp.mpf('0.05')  # ring radius R = core_fraction * cluster spacing a
as_=[mp.mpf(x) for x in ['1e-1','3e-2','1e-2','3e-3','1e-3']]
target_u=mp.mpf('1e-6')

def ring_axis(z,z0,R,Gamma,sgn):
    return sgn*Gamma*R**2/(2*(R**2+(z-z0)**2)**mp.mpf('1.5'))

def pair_unit_I(a):
    R=core_fraction*a
    Gamma=I0/(pi*R**2)
    return (ring_axis(d,a/2,R,Gamma,1)+ring_axis(d,-a/2,R,Gamma,-1), R, Gamma)

def quartet_unit_I(a):
    R=core_fraction*a
    Gamma=I0/(pi*R**2)
    st=[(-3*a/2,1),(-a/2,-1),(a/2,-1),(3*a/2,1)]
    return (sum(ring_axis(d,z0,R,Gamma,s) for z0,s in st), R, Gamma)

rows=[]
pair_u=[]; quart_u=[]; pair_Ireq=[]; quart_Ireq=[]
for a in as_:
    up,R,G=pair_unit_I(a)
    uq,_,_=quartet_unit_I(a)
    # Point-impulse asymptotic coefficients for fixed I0 and d=1.
    pred_p=3*I0*a/(2*pi*d**4)
    pred_q=12*I0*a**2/(pi*d**5)
    ratio_p=up/pred_p
    ratio_q=uq/pred_q
    pair_u.append(abs(up)); quart_u.append(abs(uq))
    Ireq_p=I0*target_u/abs(up); Ireq_q=I0*target_u/abs(uq)
    pair_Ireq.append(Ireq_p); quart_Ireq.append(Ireq_q)
    rows.append({
        'a_over_d':mp.nstr(a/d,20),
        'R_over_d':mp.nstr(R/d,20),
        'pair_u_fixed_I':mp.nstr(up,30),
        'quartet_u_fixed_I':mp.nstr(uq,30),
        'pair_ratio_to_Ia_over_d4':mp.nstr(ratio_p,30),
        'quartet_ratio_to_Ia2_over_d5':mp.nstr(ratio_q,30),
        'pair_required_I_for_target_u':mp.nstr(Ireq_p,30),
        'quartet_required_I_for_target_u':mp.nstr(Ireq_q,30),
    })

# Observed log slopes on the two smallest compactness values.
def slope(vals,xs):
    return mp.log(vals[-1]/vals[-2])/mp.log(xs[-1]/xs[-2])
sp_up=slope(pair_u,as_); sp_uq=slope(quart_u,as_)
sp_Ip=slope(pair_Ireq,as_); sp_Iq=slope(quart_Ireq,as_)
if abs(sp_up-1)>mp.mpf('0.01'): raise AssertionError(('pair influence compactness slope',sp_up))
if abs(sp_uq-2)>mp.mpf('0.02'): raise AssertionError(('quartet influence compactness slope',sp_uq))
if abs(sp_Ip+1)>mp.mpf('0.01'): raise AssertionError(('pair collateral slope',sp_Ip))
if abs(sp_Iq+2)>mp.mpf('0.02'): raise AssertionError(('quartet collateral slope',sp_Iq))
# Smallest-a asymptotics should already be close to the universal point-impulse moments.
last=rows[-1]
if abs(mp.mpf(last['pair_ratio_to_Ia_over_d4'])-1)>mp.mpf('0.01'):
    raise AssertionError(('pair asymptotic normalization',last))
if abs(mp.mpf(last['quartet_ratio_to_Ia2_over_d5'])-1)>mp.mpf('0.02'):
    raise AssertionError(('quartet asymptotic normalization',last))

print(json.dumps({
    'precision_bits_requested':BITS,
    'mpmath_dps':mp.mp.dps,
    'status':'PASS',
    'pair_influence_slope_vs_compactness':mp.nstr(sp_up,35),
    'quartet_influence_slope_vs_compactness':mp.nstr(sp_uq,35),
    'pair_required_impulse_slope':mp.nstr(sp_Ip,35),
    'quartet_required_impulse_slope':mp.nstr(sp_Iq,35),
    'interpretation':'For fixed component impulse, cancelling the total impulse suppresses remote influence linearly in donor compactness L/d; cancelling the next impulse moment suppresses it quadratically. Holding remote influence fixed therefore requires collateral growing as (d/L)^n at cancellation depth n.',
    'rows':rows,
},indent=2))
