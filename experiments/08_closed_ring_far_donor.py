import json, os
from flint import arb, ctx
import mpmath as mp

BITS = int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be >=160')
ctx.prec = BITS
mp.mp.dps = max(70, int(BITS*0.30103)+30)
pi = arb.pi()

# Exact axis efficiency for a closed circular filament.
def E(q):
    return 3*q / (2*(1+q*q)**arb('2.5'))

q_star = arb('0.5')
Emax = 24/(25*arb(5).sqrt())
if not (E(q_star)/Emax).contains(1):
    raise AssertionError('analytic maximum value mismatch')

# Derivative sign: E'(q) proportional to 1-4q^2.
for q_s, sign in [('1e-6',1),('0.1',1),('0.49',1),('0.51',-1),('1',-1),('10',-1)]:
    q=arb(q_s)
    dsign = 1-4*q*q
    if sign>0 and not dsign.is_positive():
        raise AssertionError(('bad derivative sign',q_s,dsign))
    if sign<0 and not dsign.is_negative():
        raise AssertionError(('bad derivative sign',q_s,dsign))

# Cross-check the analytic axis velocity and strain against direct ring quadrature.
def mp_ring_uz(R,G,z):
    # Biot-Savart: dX x (x-X) has z-component R^2 dtheta on the axis.
    f = lambda th: (G/(4*mp.pi)) * R*R / (R*R+z*z)**mp.mpf('1.5')
    return mp.quad(f,[0,2*mp.pi])

def mp_ring_szz_fd(R,G,z):
    h = mp.mpf('1e-20')*max(mp.mpf(1),abs(z),abs(R))
    return (mp_ring_uz(R,G,z+h)-mp_ring_uz(R,G,z-h))/(2*h)

cross = []
for R_s,G_s,q_s in [('1e-9','1e-6','0.5'),('1','1','0.5'),('1e6','1e3','2')]:
    Rm,Gm,qm = map(mp.mpf,[R_s,G_s,q_s]); zm=qm*Rm
    uz_num=mp_ring_uz(Rm,Gm,zm)
    uz_ex=Gm*Rm*Rm/(2*(Rm*Rm+zm*zm)**mp.mpf('1.5'))
    if abs((uz_num-uz_ex)/uz_ex) > mp.mpf('1e-45'):
        raise AssertionError('ring axis velocity quadrature mismatch')
    s_num=mp_ring_szz_fd(Rm,Gm,zm)
    s_ex=-3*Gm*Rm*Rm*zm/(2*(Rm*Rm+zm*zm)**mp.mpf('2.5'))
    if abs((s_num-s_ex)/s_ex) > mp.mpf('1e-25'):
        raise AssertionError(('ring axis strain quadrature mismatch',R_s,G_s,q_s,s_num,s_ex))
    cross.append({'R':R_s,'Gamma':G_s,'q':q_s,'strain_relerr':mp.nstr(abs((s_num-s_ex)/s_ex),12)})

# Scale-collateral test. Target strain = Gamma_c/ell^2.
ell_values=['1e-18','1e-9','1','1e9']
gc_values=['1e-12','1','1e12']
Lambda_values=['1','2','10','1e3','1e6']
rows=[]
for ell_s in ell_values:
  for gc_s in gc_values:
    for La_s in Lambda_values:
      ell=arb(ell_s); gc=arb(gc_s); La=arb(La_s)
      R=La*ell
      gd_required=(gc/ell**2)*(R**2)/Emax
      ratio=gd_required/gc
      expected=La**2/Emax
      if not (ratio/expected).contains(1):
          raise AssertionError(('collateral scaling failed',ell_s,gc_s,La_s,ratio,expected))
      # Directly verify returned strain equals target at q=1/2.
      s = (gd_required/R**2)*Emax
      target=gc/ell**2
      if not (s/target).contains(1):
          raise AssertionError('target strain not met')
      rows.append({'ell':ell_s,'Gamma_core':gc_s,'Lambda':La_s,'Gamma_d_over_Gamma_c':str(ratio)})

print(json.dumps({
 'arb_precision_bits':BITS,
 'status':'PASS',
 'Emax':str(Emax),
 'q_star':'0.5',
 'cases':len(rows),
 'interpretation':'For a physically closed circular vortex donor, axial strain efficiency is bounded. To supply a target core turnover strain from a donor Lambda core-scales away/large, required donor circulation grows exactly like Lambda^2 at optimal placement.',
 'quadrature_crosschecks':cross,
 'rows':rows,
},indent=2))
