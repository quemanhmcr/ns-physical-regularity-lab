import json, os, math
import mpmath as mp

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
mp.mp.dps=int(BITS*math.log10(2))+35
pi=mp.pi
alpha=mp.mpf('0.15')  # regularized core radius / ring radius; held fixed under dilation

# Base signed closed-ring stacks. Positions are in units of the ring radius.
a0=mp.mpf('3')
stacks=[
    ('single',0,[(mp.mpf('0'),1)]),
    ('pair',1,[(a0/2,1),(-a0/2,-1)]),
    ('quartet',2,[(-3*a0/2,1),(-a0/2,-1),(a0/2,-1),(3*a0/2,1)]),
]
Ls=[mp.mpf(x) for x in ['1','0.3','0.1','0.03','0.01']]

def qint(f):
    return mp.quad(f,[0,pi/2,pi,3*pi/2,2*pi])

def J(R,delta,h):
    # Double line integral reduced by rotational symmetry:
    # int int dX.dX'/sqrt(|X-X'|^2+delta^2)
    f=lambda th: mp.cos(th)/mp.sqrt(2*R**2*(1-mp.cos(th))+h**2+delta**2)
    return 2*pi*R**2*qint(f)

def energy_stack(L,Gamma,base_stack):
    R=L; delta=alpha*L
    z=[L*p for p,s in base_stack]; signs=[mp.mpf(s) for p,s in base_stack]
    cache={}
    total=mp.mpf('0')
    for i in range(len(z)):
        for j in range(len(z)):
            h=abs(z[i]-z[j])
            key=mp.nstr(h/L,40)  # normalized separation; stable cache key
            if key not in cache:
                cache[key]=J(R,delta,h)
            total += signs[i]*signs[j]*cache[key]
    return Gamma**2*total/(8*pi)

def mu(base_stack,n):
    return sum(mp.mpf(s)*p**n for p,s in base_stack)

rows=[]
for name,n,stack in stacks:
    # Verify lower signed impulse moments vanish exactly when they are supposed to.
    lower=[mu(stack,k) for k in range(n)]
    if any(v!=0 for v in lower):
        raise AssertionError((name,'lower moments not cancelled',lower))
    mun=mu(stack,n)
    if mun==0: raise AssertionError((name,'surviving moment vanished'))
    Es=[]; invariants=[]; gammas=[]
    for L in Ls:
        # Ring impulse is pi Gamma L^2.  Fix the first surviving spatial impulse moment to magnitude 1:
        # M_n = pi Gamma L^(n+2) mu_n = 1.
        Gamma=1/(pi*abs(mun)*L**(n+2))
        E=energy_stack(L,Gamma,stack)
        if E<=0: raise AssertionError((name,'regularized kinetic energy nonpositive',L,E))
        inv=E*L**(2*n+3)
        Es.append(E); invariants.append(inv); gammas.append(Gamma)
    # Exact dilation predicts E ~ L^{-(2n+3)}.
    obs=mp.log(Es[-1]/Es[-2])/mp.log(Ls[-1]/Ls[-2])
    expected=-(2*n+3)
    if abs(obs-expected)>mp.mpf('1e-18'):
        raise AssertionError((name,'packing energy exponent mismatch',obs,expected))
    inv0=invariants[0]
    max_rel=max(abs(v-inv0)/abs(inv0) for v in invariants)
    if max_rel>mp.mpf('1e-25'):
        raise AssertionError((name,'dilation invariant drift',max_rel))
    rows.append({
        'name':name,
        'cancellation_depth_n':n,
        'surviving_signed_moment_mu_n':mp.nstr(mun,30),
        'predicted_energy_exponent':expected,
        'observed_energy_exponent':mp.nstr(obs,40),
        'packing_invariant_E_times_L_power':mp.nstr(inv0,40),
        'max_relative_invariant_drift':mp.nstr(max_rel,20),
        'Gamma_values_for_fixed_Mn':[mp.nstr(g,24) for g in gammas],
        'E_values':[mp.nstr(e,24) for e in Es],
    })

print(json.dumps({
    'precision_bits_requested':BITS,
    'mpmath_dps':mp.mp.dps,
    'status':'PASS',
    'regularization_core_ratio':mp.nstr(alpha,20),
    'interpretation':'For geometrically similar regularized closed-vortex stacks with the first surviving impulse moment held fixed, kinetic energy scales as L^-(2n+3): L^-3 for impulse, L^-5 after impulse cancellation, and L^-7 after one further moment cancellation. This is a canonical impulse-packing energy ladder, not yet a universal lower bound for arbitrary Navier-Stokes vorticity.',
    'rows':rows,
},indent=2))
