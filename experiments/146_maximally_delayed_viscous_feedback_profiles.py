import json, os
from fractions import Fraction as Q
from flint import arb,ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160:raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def M(l,k):
 # C_l[r^(2k)] at L=1.
 return -Q(l+1,(2*k+2)*(2*k+2*l+3))
def fac(l,k,j):
 # coefficient of D_l^j r^(2k), j<=k.
 q=Q(1)
 for s in range(j):
  p=k-s;q*=2*p*(2*p+2*l+1)
 return q
def solve_delayed(l,n):
 # c_n=1. Endpoint constraints E_j=0, j=0..n-2 are triangular downward.
 c=[Q(0)]*(n+1);c[n]=Q(1)
 for j in range(n-2,-1,-1):
  # E_j = -sum_{k>j} c_k fac(k,j)=0. Solve c_{j+1}; higher terms already known.
  denom=fac(l,j+1,j);rest=sum((c[k]*fac(l,k,j) for k in range(j+2,n+1)),Q(0));c[j+1]=-rest/denom
 # screened moment determines c0.
 c[0]=-sum((c[k]*M(l,k) for k in range(1,n+1)),Q(0))/M(l,0)
 return c
def Citer(l,c,j):
 # C_l[D^j a]. j=0 uses screened moment; j>=1 endpoint law.
 if j==0:return sum((c[k]*M(l,k) for k in range(len(c))),Q(0))
 # (l+1)(D^(j-1)a(0)-D^(j-1)a(1)) = -(l+1) sum_{k>=j} c_k fac(k,j-1)
 return -Q(l+1)*sum((c[k]*fac(l,k,j-1) for k in range(j,len(c))),Q(0))
def radial_norm(l,c):
 # angular factor cancels; integral r^(2l+2) a^2 dr.
 return sum((c[i]*c[j]*Q(1,2*l+2+2*i+2*j+1) for i in range(len(c)) for j in range(len(c))),Q(0))
def radial_dirichlet(l,c):
 # integral r^(2l+2) |a'|^2 dr.
 return sum((c[i]*c[j]*(2*i)*(2*j)*Q(1,2*l+2+2*i+2*j-1) for i in range(1,len(c)) for j in range(1,len(c))),Q(0))
def A(q):return arb(q.numerator)/q.denominator
rows=[]
for l in (2,4,6,8):
 for n in range(1,9):
  c=solve_delayed(l,n);vals=[Citer(l,c,j) for j in range(n+1)]
  for j in range(n):
   if vals[j]!=0:raise AssertionError(('delay constraint',l,n,j,vals[j]))
  if vals[n]==0:raise AssertionError(('first reveal vanished',l,n))
  Z=radial_norm(l,c);D=radial_dirichlet(l,c)
  if Z<=0 or D<=0:raise AssertionError(('physical forms',l,n,Z,D))
  ray=D/Z;reveal=vals[n]*vals[n]/Z
  rows.append({'l':l,'maximal_hidden_viscous_orders':n-1,'polynomial_radial_degree':2*n,'coefficients_c0_to_cn':[str(A(x)) for x in c],'first_nonzero_feedback_derivative_order':n,'first_reveal_C_l_Dn_a':str(A(vals[n])),'radial_enstrophy_weight_integral':str(A(Z)),'radial_Dirichlet_weight_integral':str(A(D)),'scale_invariant_Dirichlet_over_enstrophy':str(A(ray)),'scale_invariant_first_reveal_square_over_enstrophy':str(A(reveal))})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For every l and every finite hiding depth n-1, there is a unique (up to scale) even polynomial radial profile of degree 2n satisfying C_l[D_l^j a]=0 for j=0,...,n-1 and C_l[D_l^n a] nonzero.  The endpoint constraints are triangular and the screened moment fixes the constant term.  Thus no fixed finite number of viscous derivatives detects every hidden radial profile.  The physically meaningful burden is instead measured by the scale-invariant radial Dirichlet/enstrophy quotient and by the first revealed feedback amplitude normalized by radial enstrophy; these are reported without treating polynomial order itself as a resource.','rows':rows},indent=2,allow_nan=False))
