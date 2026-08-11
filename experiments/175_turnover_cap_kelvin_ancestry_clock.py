import json, os
from fractions import Fraction as F
from math import comb
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import degree6_hodge_servo_core as C
z=C.z;one=C.o;pi=arb.pi();rt3=arb(3).sqrt();K=arb(4)*pi/(3*rt3)
X=({(1,0,0):one},{(0,1,0):one},{(0,0,1):one})
B=((arb(2),z,z),(z,-one,z),(z,z,-one))

def ppow(P,n):
 q={(0,0,0):one}
 for _ in range(n):q=C.pmul(q,P)
 return q
r2=C.padd(C.padd(C.pmul(X[0],X[0]),C.pmul(X[1],X[1])),C.pmul(X[2],X[2]))
omega0=[]
for i in range(3):
 q={}
 for j in range(3):q=C.padd(q,C.pscale(B[i][j],X[j]))
 omega0.append(q)
omega0=tuple(omega0);u0=C.vscale(-arb(1)/3,C.cross(X,omega0))
rows=[]
for p in (2,3,4,8,16):
 # chi=(1-r^2)^p at L=1.
 chi={}
 for k in range(p+1): chi=C.padd(chi,C.pscale(((-1)**k)*comb(p,k),ppow(r2,k)))
 u=tuple(C.pmul(chi,q) for q in u0)
 lapu=tuple(C.plap(q) for q in u)
 # Structural radial operator Dchi=chi''+6 chi'/r on a degree-2 harmonic velocity.
 Dchi={}
 for k in range(1,p+1):
  ck=arb(((-1)**k)*comb(p,k));lam=arb(2*k*(2*k+5));Dchi=C.padd(Dchi,C.pscale(ck*lam,ppow(r2,k-1)))
 pred=tuple(C.pmul(Dchi,q) for q in u0)
 err=C.savg(C.vdot(C.vadd(lapu,C.vscale(-1,pred)),C.vadd(lapu,C.vscale(-1,pred))))
 if not err.contains(0):raise AssertionError(('vector Laplacian radial identity',p,err))
 speak=(arb(3)/arb(2*p+3)).sqrt(); s2=speak*speak
 chipeak=(1-s2)**p
 D_over_chi=-(arb(8)*p+18+arb(9)/p)
 # Independent direct formula at peak.
 chip=-2*p*speak*(1-s2)**(p-1)
 chipp=-2*p*(1-s2)**(p-1)+4*p*(p-1)*s2*(1-s2)**(p-2)
 direct=(chipp+6*chip/speak)/chipeak
 if not (direct/D_over_chi).contains(1):raise AssertionError(('peak Kelvin coefficient',p,direct,D_over_chi))
 # Material cap is fixed as a set under u_chi. Kelvin gives Gamma_dot/Gamma=nu/L^2 Dchi/chi.
 rows.append({
   'taper_power_p':p,'peak_cap_radius_over_L':str(speak),'chi_at_peak':str(chipeak),
   'fractional_Kelvin_circulation_rate_coefficient_times_nu_over_L2':str(D_over_chi),
   'instantaneous_Kelvin_ancestry_decay_clock_coefficient_times_L2_over_nu':str(-1/D_over_chi),
   'vector_Laplacian_identity_error_square':str(err),
   'direct_peak_Dchi_over_chi':str(direct),
 })
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'For the finite-energy localized null catalyst u=chi(r/L)u0, every centered spherical cap remains material under the catalyst velocity because u is purely azimuthal around the x-axis.  Kelvin therefore measures ancestry mutation directly on the cap boundary. '
  'The uncut degree-two velocity u0 is componentwise harmonic.  Exact polynomial certification gives Delta(chi u0)=[chi_rr+6 chi_r/r]u0.  Hence the material cap circulation satisfies Gamma_dot_visc/Gamma=(nu/L^2)[chi_ss+6 chi_s/s]/chi. '
  'For chi=(1-s^2)^p the radial cap flux peaks at s^2=3/(2p+3), and the exact fractional Kelvin rate there is -[8p+18+9/p] nu/L^2.  Thus the turnover collar is not merely an enstrophy observer: it is precisely where the reusable Euler circulation ancestry begins to be viscously destroyed on an order-L^2/nu clock. '
  'Combined with module173, this clock is still asymptotically longer than the time remaining in the critical branch, so viscosity cannot by itself create the required divergent cap flux near T.  The remaining positive flux growth must be supplied by source-relative material recruitment/external deformation.'),
 'rows':rows
},indent=2,allow_nan=False))
