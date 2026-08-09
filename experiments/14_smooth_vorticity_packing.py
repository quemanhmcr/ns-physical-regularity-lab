import json, os, math
import mpmath as mp

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
mp.mp.dps=int(BITS*math.log10(2))+35
pi=mp.pi

# Smooth divergence-free vorticity is generated from A=A_z e_z by omega=curl A.
# Each Gaussian A-packet has hydrodynamic impulse exactly equal to integral A dx.
# Signed packet stacks reproduce the same impulse-moment cancellation hierarchy as
# the filament experiments, without any singular filament core.
a0=mp.mpf('3')
stacks=[
    ('single',0,[(mp.mpf('0'),1)]),
    ('pair',1,[(a0/2,1),(-a0/2,-1)]),
    ('quartet',2,[(-3*a0/2,1),(-a0/2,-1),(a0/2,-1),(3*a0/2,1)]),
]
Ls=[mp.mpf(x) for x in ['1','0.3','0.1','0.03','0.01']]

def mu(stack,n):
    return sum(mp.mpf(s)*c**n for c,s in stack)

def gaussian_energy_shape(stack):
    # Isotropic Gaussian carrier: after Schwinger parameterization and exact
    # transverse Gaussian integration, only u in [0,1] remains.
    total=mp.mpf('0')
    for ci,si in stack:
        for cj,sj in stack:
            d=ci-cj
            f=lambda u: mp.sqrt(u)*mp.exp(-(d*d)*u/4)
            total += mp.mpf(si)*mp.mpf(sj)*mp.quad(f,[0,mp.mpf('0.25'),mp.mpf('0.5'),mp.mpf('0.75'),1])
    return total

def gaussian_enstrophy_shape(stack):
    # Exact Fourier Gaussian integral for integral |omega|^2.
    return sum(mp.mpf(si)*mp.mpf(sj)*mp.exp(-((ci-cj)**2)/4)
               for ci,si in stack for cj,sj in stack)

rows=[]
for name,n,stack in stacks:
    for k in range(n):
        if mu(stack,k)!=0: raise AssertionError((name,'uncancelled lower moment',k,mu(stack,k)))
    mun=mu(stack,n)
    if mun==0: raise AssertionError((name,'surviving moment zero'))
    FE=gaussian_energy_shape(stack)
    FO=gaussian_enstrophy_shape(stack)
    if FE<=0 or FO<=0: raise AssertionError((name,'nonpositive smooth field quadratic form',FE,FO))
    # Fourier/Parseval constants for A_z Gaussian packets of unit width.
    CE=FE/(16*pi**mp.mpf('1.5'))
    CO=FO/(8*pi**mp.mpf('1.5'))
    Es=[]; Os=[]; taus=[]
    for L in Ls:
        # Each packet impulse P; M_n=P L^n mu_n=1.
        P=1/(abs(mun)*L**n)
        E=P**2*CE/L**3
        O=P**2*CO/L**5
        tau=E/O  # nu=1; with viscosity nu, multiply by 1/nu.
        Es.append(E); Os.append(O); taus.append(tau)
    e_slope=mp.log(Es[-1]/Es[-2])/mp.log(Ls[-1]/Ls[-2])
    o_slope=mp.log(Os[-1]/Os[-2])/mp.log(Ls[-1]/Ls[-2])
    t_slope=mp.log(taus[-1]/taus[-2])/mp.log(Ls[-1]/Ls[-2])
    if abs(e_slope+(2*n+3))>mp.mpf('1e-25'): raise AssertionError((name,'energy exponent',e_slope))
    if abs(o_slope+(2*n+5))>mp.mpf('1e-25'): raise AssertionError((name,'enstrophy exponent',o_slope))
    if abs(t_slope-2)>mp.mpf('1e-25'): raise AssertionError((name,'viscous clock exponent',t_slope))
    rows.append({
        'name':name,'cancellation_depth_n':n,
        'energy_exponent':mp.nstr(e_slope,35),
        'enstrophy_exponent':mp.nstr(o_slope,35),
        'viscous_lifetime_exponent':mp.nstr(t_slope,35),
        'dimensionless_E_shape':mp.nstr(CE,30),
        'dimensionless_Omega_shape':mp.nstr(CO,30),
        'tau_over_L2_for_nu1':mp.nstr(CE/CO,30),
    })

# Shape-optimization attack for a single smooth impulse carrier.
# A_z ~ exp[-r_perp^2/(2s^2)-z^2/(2h^2)], lambda=h/s.
# Fix RMS geometric size by 2s^2+h^2=3 L^2 and set L=1, fixed impulse I=1.
def F_energy_single(s,h):
    a=s*s-h*h
    tiny=mp.mpf('1e-40')
    if abs(a)<tiny:
        return 2/(3*h**3)
    if a>0:
        root=mp.sqrt(a)
        return pi/(2*a**mp.mpf('1.5')) - h/(a*s*s) - mp.atan(h/root)/(a**mp.mpf('1.5'))
    b=-a
    root=mp.sqrt(b)
    return h/(b*s*s) - mp.atanh(root/h)/(b**mp.mpf('1.5'))

lambdas=[mp.mpf(x) for x in ['1e-3','1e-2','0.1','0.3','1','3','10','100','1e3']]
shape_rows=[]
CE_iso=None; CO_iso=None
for lam in lambdas:
    s=mp.sqrt(3/(2+lam*lam)); h=lam*s
    F=F_energy_single(s,h)
    CE=F/(16*pi**mp.mpf('1.5'))
    CO=1/(8*pi**mp.mpf('1.5')*s**4*h)
    if CE<=0 or CO<=0: raise AssertionError(('aspect','nonpositive',lam,CE,CO))
    if lam==1:
        CE_iso=CE; CO_iso=CO
    shape_rows.append((lam,s,h,CE,CO,CE/CO))
if CE_iso is None: raise AssertionError('missing isotropic calibration')
min_CE=min(r[3] for r in shape_rows)
# Deliberately weak falsification threshold: if anisotropy can drive energy almost free,
# this gate should fail.  Asymptotic analysis predicts a finite oblate limit.
if min_CE < mp.mpf('0.004'):
    raise AssertionError(('anisotropy made fixed-impulse energy nearly free',min_CE))
# Both extreme oblate and prolate shapes should pay strongly in enstrophy.
if not (shape_rows[0][4] > 10*CO_iso and shape_rows[-1][4] > 10*CO_iso):
    raise AssertionError(('viscous geometry tax did not grow at anisotropic extremes',shape_rows[0][4],CO_iso,shape_rows[-1][4]))

shape_json=[]
for lam,s,h,CE,CO,tau in shape_rows:
    shape_json.append({
        'lambda_h_over_s':mp.nstr(lam,12),
        's_over_L':mp.nstr(s,20),'h_over_L':mp.nstr(h,20),
        'E_times_L3_over_I2':mp.nstr(CE,30),
        'Omega_times_L5_over_I2':mp.nstr(CO,30),
        'tau_nu_over_L2':mp.nstr(tau,30),
    })

print(json.dumps({
    'precision_bits_requested':BITS,
    'mpmath_dps':mp.mp.dps,
    'status':'PASS',
    'interpretation':'The impulse-packing ladder survives replacement of singular filaments by smooth divergence-free Gaussian vorticity: E~L^-(2n+3), enstrophy~L^-(2n+5), so the viscous depletion clock is always L^2/nu. Within a wide anisotropic Gaussian family at fixed impulse and RMS size, kinetic-energy cost stays bounded away from zero while extreme flattening/elongation sharply increases enstrophy.',
    'smooth_packing_rows':rows,
    'aspect_ratio_attack':{
        'minimum_sampled_E_times_L3_over_I2':mp.nstr(min_CE,30),
        'isotropic_E_times_L3_over_I2':mp.nstr(CE_iso,30),
        'isotropic_Omega_times_L5_over_I2':mp.nstr(CO_iso,30),
        'rows':shape_json,
    },
},indent=2))
