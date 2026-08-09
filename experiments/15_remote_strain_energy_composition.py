import json, os, math
import mpmath as mp

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
mp.mp.dps=int(BITS*math.log10(2))+35
pi=mp.pi
alpha=mp.mpf('0.15')

a0=mp.mpf('3')
stacks=[
    ('single',0,[(mp.mpf('0'),1)]),
    ('pair',1,[(a0/2,1),(-a0/2,-1)]),
    ('quartet',2,[(-3*a0/2,1),(-a0/2,-1),(a0/2,-1),(3*a0/2,1)]),
]

def qint(f):
    return mp.quad(f,[0,pi/2,pi,3*pi/2,2*pi])

def J(R,delta,h):
    f=lambda th: mp.cos(th)/mp.sqrt(2*R**2*(1-mp.cos(th))+h**2+delta**2)
    return 2*pi*R**2*qint(f)

def energy_stack(L,Gamma,stack):
    R=L; delta=alpha*L
    z=[L*p for p,s in stack]; signs=[mp.mpf(s) for p,s in stack]
    cache={}; total=mp.mpf('0')
    for i in range(len(z)):
        for j in range(len(z)):
            h=abs(z[i]-z[j]); key=mp.nstr(h/L,40)
            if key not in cache: cache[key]=J(R,delta,h)
            total += signs[i]*signs[j]*cache[key]
    return Gamma**2*total/(8*pi)

def ring_axis_strain(d,z0,R,Gamma,sgn):
    q=d-z0
    return -3*mp.mpf(sgn)*Gamma*R**2*q/(2*(R**2+q**2)**mp.mpf('2.5'))

def stack_strain(d,L,Gamma,stack):
    return sum(ring_axis_strain(d,L*c,L,Gamma,s) for c,s in stack)

def slope(vals,xs):
    return mp.log(vals[-1]/vals[-2])/mp.log(xs[-1]/xs[-2])

target_s=mp.mpf('1')
etas=[mp.mpf(x) for x in ['0.1','0.05','0.03','0.02','0.01','0.005']]
compact_rows=[]
for name,n,stack in stacks:
    qs=[]; gammas=[]
    for eta in etas:
        d=mp.mpf('1'); L=eta*d
        su=stack_strain(d,L,1,stack)
        if su==0: raise AssertionError((name,'zero unit strain',eta))
        Gamma=target_s/abs(su)
        E=energy_stack(L,Gamma,stack)
        Q=E/(target_s**2*d**5)
        qs.append(Q); gammas.append(Gamma)
    obs=slope(qs,etas); expected=-(2*n+3)
    tol=[mp.mpf('0.03'),mp.mpf('0.08'),mp.mpf('0.16')][n]
    if abs(obs-expected)>tol:
        raise AssertionError((name,'compactness composition exponent',obs,expected))
    compact_rows.append({
        'name':name,'cancellation_depth_n':n,
        'predicted_compactness_exponent':expected,
        'observed_compactness_exponent':mp.nstr(obs,35),
        'eta_values':[mp.nstr(x,16) for x in etas],
        'E_over_s2d5':[mp.nstr(x,25) for x in qs],
        'required_Gamma':[mp.nstr(x,25) for x in gammas],
    })

# With compactness eta=L/d held fixed, exact NS/Biot-Savart scaling predicts
# E_required ~ s^2 d^5 for every cancellation depth.
eta0=mp.mpf('0.02')
ds=[mp.mpf(x) for x in ['1e-6','1e-3','1','1e3','1e6']]
distance_rows=[]
for name,n,stack in stacks:
    Es=[]; Gs=[]; invariants=[]
    for d in ds:
        L=eta0*d
        su=stack_strain(d,L,1,stack)
        Gamma=target_s/abs(su)
        E=energy_stack(L,Gamma,stack)
        Es.append(E); Gs.append(Gamma); invariants.append(E/(target_s**2*d**5))
    eslope=slope(Es,ds); gslope=slope(Gs,ds)
    if abs(eslope-5)>mp.mpf('1e-24'): raise AssertionError((name,'distance energy exponent',eslope))
    if abs(gslope-2)>mp.mpf('1e-24'): raise AssertionError((name,'distance circulation exponent',gslope))
    inv0=invariants[0]
    drift=max(abs(x-inv0)/abs(inv0) for x in invariants)
    if drift>mp.mpf('1e-24'): raise AssertionError((name,'distance invariant drift',drift))
    distance_rows.append({
        'name':name,'cancellation_depth_n':n,
        'energy_distance_exponent':mp.nstr(eslope,35),
        'required_Gamma_distance_exponent':mp.nstr(gslope,35),
        'E_over_s2d5_at_eta_fixed':mp.nstr(inv0,35),
        'max_relative_invariant_drift':mp.nstr(drift,20),
    })

print(json.dumps({
    'precision_bits_requested':BITS,'mpmath_dps':mp.mp.dps,'status':'PASS',
    'interpretation':'Exact closed-vortex stacks confirm the composed remote-collateral law. After eliminating the surviving impulse moment, the energy required to impose a fixed remote strain scales universally as s^2 d^5 at fixed source compactness, independent of cancellation depth. Cancelling n lower impulse moments survives only as the compactness penalty (d/L)^(2n+3). Required circulation scales as s d^2 at fixed L/d for every depth.',
    'compactness_attack':compact_rows,
    'distance_scaling':distance_rows,
},indent=2))
