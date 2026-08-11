import json, os, math
from fractions import Fraction as Q
from flint import arb,ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160:raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def M(l,k):return -Q(l+1,(2*k+2)*(2*k+2*l+3))
def fac(l,k,j):
 q=Q(1)
 for s in range(j):
  p=k-s;q*=2*p*(2*p+2*l+1)
 return q
def coeffs(l,n):
 c=[Q(0)]*(n+1);c[n]=Q(1)
 for j in range(n-2,-1,-1):c[j+1]=-sum((c[k]*fac(l,k,j) for k in range(j+2,n+1)),Q(0))/fac(l,j+1,j)
 c[0]=-sum((c[k]*M(l,k) for k in range(1,n+1)),Q(0))/M(l,0);return c
def Citer(l,c,j):
 if j==0:return sum((c[k]*M(l,k) for k in range(len(c))),Q(0))
 return -Q(l+1)*sum((c[k]*fac(l,k,j-1) for k in range(j,len(c))),Q(0))
def Z(l,c):return sum((c[i]*c[j]*Q(1,2*l+2+2*i+2*j+1) for i in range(len(c)) for j in range(len(c))),Q(0))
def A(q):return arb(q.numerator)/q.denominator
rows=[]
for l in (2,4,6,8):
 for n in range(1,17):
  c=coeffs(l,n);z=Z(l,c);v=Citer(l,c,n)
  radial_amp=abs(A(v))/A(z).sqrt()
  hodge_amp=abs(A(v))/(arb(l+1)*A(z)).sqrt()
  radial_theta=(arb(math.factorial(n))/radial_amp)**(arb(1)/n)
  theta=(arb(math.factorial(n))/hodge_amp)**(arb(1)/n)
  # exact leading coefficient check from the top monomial c_n=1.
  closed=M(l,0)*fac(l,n,n)
  if v!=closed:raise AssertionError(('top reveal closed form',l,n,v,closed))
  rows.append({'l':l,'hiding_order_n':n,'radial_coordinate_first_reveal_amplitude_abs':str(radial_amp),'radial_coordinate_reveal_clock':str(radial_theta),'physical_Hodge_energy_first_reveal_amplitude_abs':str(hodge_amp),'dimensionless_exact_Hodge_energy_reveal_clock_Theta':str(theta),'n_times_Theta':str(theta*n),'four_n_times_Theta':str(theta*4*n),'exact_C_l_Dn_a':str(A(v))})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'Write s=r/L and omega=Omega b(s) X_l(sphere), so tau=nu t/L^2 is the exact heat time.  For the maximally delayed degree-2n radial polynomial, all screened feedback heat coefficients below order n vanish and D_l^(n+1)b=0, hence the source-scaled feedback is exactly tau^n C_l[D_l^n b]/n!.  The physical Hodge-energy normalization satisfies 2 E_h/(L^2 Zeta_0)=|F_l|^2/((l+1) Z), where Z is the reported radial enstrophy factor.  Therefore Theta=[n! sqrt((l+1)Z)/|C_l[D_l^n b]|]^(1/n) is the exact dimensionless time at which the generated lower harmonic kinetic occupancy reaches 2 E_h=L^2 Zeta_0.  A radial-coordinate clock without the angular Hodge factor is also reported only as an observer diagnostic.  The first reveal coefficient is independently checked from the top radial monomial.','rows':rows},indent=2,allow_nan=False))
