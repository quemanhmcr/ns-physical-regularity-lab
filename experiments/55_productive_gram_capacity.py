import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
sqrt3=arb(3).sqrt(); cap=arb(4)/27; one=arb(1)

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x > 1-t and x < 1+t): raise AssertionError((label,x))

# Exact global equality geometry.
r=1/sqrt3; s=(arb(2)/3).sqrt()
a=(s,arb(0),r); n=(arb(0),arb(0),one); b=(arb(0),-s,-r)
alpha=dot(a,n); beta=dot(b,n); gamma=dot(a,b); T=dot(a,cross(n,b))
F=1+2*alpha*beta*gamma-alpha*alpha-beta*beta-gamma*gamma
P=(-alpha*beta)*F
Kba=T*alpha; Kab=-T*beta
certify_one(F/(T*T),('Gram determinant equals T2',))
certify_one(P/cap,('global capacity equality',))
certify_one(Kba/(arb(2)/(3*sqrt3)),('edge equality Kba',)); certify_one(Kab/(arb(2)/(3*sqrt3)),('edge equality Kab',))

# Exact one-dimensional factorization underlying the global bound.
xs=['0','1e-24','0.01','0.1','0.5773502691896257645','0.9','0.999999','1']
# Use an Arb-native exact optimum point separately; decimal samples only attack the factorization away from it.
rows=[]
M=arb(2)/(3*sqrt3)
for xs_ in xs:
    x=arb(xs_)
    lhs=M-x*(1-x*x)
    rhs=(x-r)*(x-r)*(x+2*r)
    if not (lhs-rhs).contains(0): raise AssertionError(('capacity factorization',xs_,lhs,rhs))
    if not (lhs>=0): raise AssertionError(('negative capacity gap',xs_,lhs))
    rows.append({'x':xs_,'x1minusx2':str(x*(1-x*x)),'capacity_gap':str(lhs),'factorized_gap':str(rhs)})

# Broad direct triad attack: alpha=x>0, beta=-y<0, transverse angle phi.
vals=['0.01','0.1','0.3','0.5773502691896257','0.9','0.99']
phis=['0.1','0.7','1.5707963267948966']
attacks=[]
for xs_ in vals:
  x=arb(xs_); ax=(1-x*x).sqrt()
  for ys_ in vals:
    y=arb(ys_); by=(1-y*y).sqrt()
    for ps in phis:
      ph=arb(ps); cp=ph.cos(); sp=ph.sin()
      aa=(ax,arb(0),x); bb=(by*cp,-by*sp,-y)
      al=dot(aa,n); be=dot(bb,n); ga=dot(aa,bb); tt=dot(aa,cross(n,bb))
      ff=1+2*al*be*ga-al*al-be*be-ga*ga; pp=(-al*be)*ff
      if not (ff>=0 and pp>=0 and pp<=cap): raise AssertionError(('Gram capacity attack',xs_,ys_,ps,ff,pp,cap))
      if not (ff-tt*tt).contains(0): raise AssertionError(('Gram identity attack',xs_,ys_,ps,ff,tt*tt))
      attacks.append({'alpha_abs':xs_,'beta_abs':ys_,'phi':ps,'T2':str(ff),'P':str(pp),'P_over_4over27':str(pp/cap)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','factor_cases':len(rows),'triad_cases':len(attacks),
 'equality':{'alpha':str(alpha),'beta':str(beta),'gamma':str(gamma),'T':str(T),'K_b_to_a':str(Kba),'K_a_to_b':str(Kab),'P':str(P),'P_over_capacity':str(P/cap)},
 'interpretation':(
  'The intrinsic triad Gram determinant equals T^2, and the full positive pair-cycle product is (-alpha beta) det G. '
  'The exact factorization 2/(3sqrt3)-x(1-x^2)=(x-1/sqrt3)^2(x+2/sqrt3) yields the global capacity P<=4/27. '
  'The symmetric geometry alpha=-beta=1/sqrt3 with perpendicular transverse projections saturates the bound and both edge weights equal 2/(3sqrt3).'
 ),'factorization':rows,'attacks':attacks
},indent=2,allow_nan=False))
