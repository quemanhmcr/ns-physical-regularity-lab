import json, os, math
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

# Exact polynomial calibration of the thin-toroidal single-loop attack.
# On |X|,|Z|<=1 use phi(s)=(1-s^2)^k, extended by zero.  With k=8 the
# profile is C^7 at the support edge.  The theory uses a C-infinity bump;
# this high-order polynomial lets every volume integral below be evaluated
# algebraically and exhibits the same scaling.
#
# u_theta=A*(r/R)*phi(X)*phi(Z), X=(r-R)/delta, Z=z/delta.
# At C: r=R,z=0, Delta_XZ[phi(X)phi(Z)](0)=-4k, hence
# |integral_C Delta u.dl| = 8*pi*k*R*A/delta^2.
# Choose A=delta^2 (R=1): Kelvin defect is constant while enstrophy ~delta^4.

k=8
R=arb(1)

def padd(a,b):
    n=max(len(a),len(b)); out=[arb(0) for _ in range(n)]
    for i,x in enumerate(a): out[i]+=x
    for i,x in enumerate(b): out[i]+=x
    return out

def pscale(a,c): return [x*c for x in a]
def pmul(a,b):
    out=[arb(0) for _ in range(len(a)+len(b)-1)]
    for i,x in enumerate(a):
      for j,y in enumerate(b): out[i+j]+=x*y
    return out

def pder(a): return [arb(i)*a[i] for i in range(1,len(a))]
def pint_sym(a):
    s=arb(0)
    for n,c in enumerate(a):
      if n%2==0: s += c*2/arb(n+1)
    return s

# phi=(1-X^2)^k
phi=[arb(0) for _ in range(2*k+1)]
for j in range(k+1):
    phi[2*j]=arb(((-1)**j)*math.comb(k,j))
dphi=pder(phi)
phi2=pmul(phi,phi); dphi2=pmul(dphi,dphi)
Iz_phi2=pint_sym(phi2); Iz_dphi2=pint_sym(dphi2)

rows=[]
prev_omega=None
prev_delta=None
scaled_tail=[]
for ds in ['0.1','0.03','0.01','0.003','0.001','1e-4','1e-6','1e-9','1e-12','1e-18','1e-24','1e-30']:
    delta=arb(ds)
    if not (delta<R/2): raise AssertionError(('delta must remain inside torus',ds))
    A=delta*delta
    rpoly=[R,delta]
    r2=pmul(rpoly,rpoly); r3=pmul(r2,rpoly)
    # omega_r contribution after change of variables:
    # Ir = int r^3 phi(X)^2 dX * int phi'(Z)^2 dZ.
    Ir=pint_sym(pmul(r3,phi2))*Iz_dphi2
    # omega_z=(A/R)[2phi + (r/delta)phi']phi(Z).
    r_over_delta=[R/delta,arb(1)]
    bracket=padd(pscale(phi,arb(2)),pmul(r_over_delta,dphi))
    Ix_z=pint_sym(pmul(rpoly,pmul(bracket,bracket)))
    Iomega_z=delta*delta*Ix_z*Iz_phi2
    Omega=2*pi*A*A/(R*R)*(Ir+Iomega_z)
    defect=8*pi*arb(k)*R*A/(delta*delta)
    target=8*pi*arb(k)*R
    if not (defect/target).contains(1):
        raise AssertionError(('Kelvin defect failed to remain constant',ds,defect,target))
    if not (Omega>0): raise AssertionError(('nonpositive enstrophy',ds,Omega))
    scaled=Omega/(delta**4)
    if prev_omega is not None and not (Omega<prev_omega):
        raise AssertionError(('enstrophy failed to decrease as layer thinned',prev_delta,ds,prev_omega,Omega))
    prev_omega=Omega; prev_delta=ds
    if delta<=arb('1e-6'): scaled_tail.append(scaled)
    rows.append({'delta':ds,'A_delta_equals_delta2':str(A),
                 'abs_Kelvin_defect':str(defect),
                 'enstrophy':str(Omega),
                 'enstrophy_over_delta4':str(scaled),
                 'defect_ratio_to_constant':str(defect/target)})

# The exact polynomial volume integral resolves the predicted delta^4 law.
for i in range(1,len(scaled_tail)):
    if not (abs(scaled_tail[i]/scaled_tail[-1]-1)<arb('1e-8')):
        raise AssertionError(('thin-torus enstrophy/delta^4 asymptotic not stable',i,scaled_tail[i],scaled_tail[-1]))

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'profile_power_k':k,
  'cases':len(rows),
  'tail_scaled_enstrophy':str(scaled_tail[-1]),
  'interpretation':'A thin divergence-free toroidal swirl around one circle can keep the instantaneous Kelvin defect exactly constant while its volume enstrophy decays like delta^4. This kills a static single-loop defect-to-enstrophy toll. The construction does not make a finite circulation renewal free in time: the layer viscous clock also collapses like delta^2/nu, so the surviving question is a spacetime finite-thickness collar-breach law.',
  'rows':rows,
},indent=2,allow_nan=False))
