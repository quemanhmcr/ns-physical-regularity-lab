import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def norm(a): return dot(a,a).sqrt()
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x > 1-t and x < 1+t):
        raise AssertionError((label,'ratio not tightly certified around one',x))

# Exact steady NS family u=(-g x+eps y z, g y,0), with c=eps/(2g)>0.
# Work entirely in strain time N=g t and dimensionless c, so no large-parent-state subtraction occurs.
cs=['1e-6','0.01','0.1','1','10','1e6']
Ls=['0.5','2','10']
Ns=['0','1e-12','0.01','0.1','1','2','5','10','30','100','200']
rows=[]
for cs_ in cs:
  c=arb(cs_)
  for Ls_ in Ls:
    L=arb(Ls_)
    for Ns_ in Ns:
      N=arb(Ns_); x=(-N).exp(); x2=x*x
      H=c+(L-c)*x2
      Rt=(-H,arb(-2),x)  # scaled bridge x R
      R2=dot(Rt,Rt); Rn=R2.sqrt()
      xa=(arb(0),arb(1),arb(0))
      xb=(arb(0),arb(-1),-x)
      xb_n=norm(xb); xb=tuple(v/xb_n for v in xb)
      n=tuple(v/Rn for v in Rt)
      T=dot(xa,cross(n,xb))
      alpha=dot(xa,n); beta=dot(xb,n)
      Kba=T*alpha; Kab=-T*beta
      if not (Kba>0 and Kab>0):
          raise AssertionError(('positive cycle sign lost',cs_,Ls_,Ns_,T,alpha,beta,Kba,Kab))

      T_closed=-x*H/((1+x2).sqrt()*Rn)
      alpha_closed=-arb(2)/Rn
      beta_closed=(2-x2)/((1+x2).sqrt()*Rn)
      Kba_closed=2*x*H/((1+x2).sqrt()*R2)
      Kab_closed=x*H*(2-x2)/((1+x2)*R2)
      certify_one(T/T_closed,('T closed',cs_,Ls_,Ns_))
      certify_one(alpha/alpha_closed,('alpha closed',cs_,Ls_,Ns_))
      certify_one(beta/beta_closed,('beta closed',cs_,Ls_,Ns_))
      certify_one(Kba/Kba_closed,('Kba closed',cs_,Ls_,Ns_))
      certify_one(Kab/Kab_closed,('Kab closed',cs_,Ls_,Ns_))

      rows.append({
        'c_eps_over_2g':cs_,'L':Ls_,'strain_time_N':Ns_,
        'x_exp_minus_N':str(x),'H':str(H),
        'T':str(T),'alpha':str(alpha),'beta':str(beta),
        'K_b_to_a':str(Kba),'K_a_to_b':str(Kab),
        'scaled_bridge_norm_xr':str(Rn),
      })

# Asymptotic coefficient: both positive edges become [2c/(c^2+4)] e^{-N}.
# At N=200 the correction is far below 1e-20 for the tested finite geometry range.
asym=[]
N=arb(200); x=(-N).exp(); x2=x*x
for cs_ in cs:
  c=arb(cs_)
  for Ls_ in Ls:
    L=arb(Ls_); H=c+(L-c)*x2; R2=H*H+4+x2
    Kba=2*x*H/((1+x2).sqrt()*R2)
    Kab=x*H*(2-x2)/((1+x2)*R2)
    coeff=2*c/(c*c+4)
    rb=Kba/(coeff*x); ra=Kab/(coeff*x)
    certify_one(rb,('Kba asymptotic coefficient',cs_,Ls_),tol='1e-20')
    certify_one(ra,('Kab asymptotic coefficient',cs_,Ls_),tol='1e-20')
    asym.append({'c':cs_,'L':Ls_,'Kba_over_coeff_exp_minus_N':str(rb),'Kab_over_coeff_exp_minus_N':str(ra),'coeff_2c_over_c2plus4':str(coeff)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'asymptotic_cases':len(asym),
 'interpretation':(
   'The exact steady quadratic Navier-Stokes family carries a material endpoint pair whose two directed angular transactions remain strictly positive for every finite strain time. '
   'Nevertheless the normalized triple-product cell T and both edge weights decay asymptotically like exp(-N); the edge coefficient is exactly 2c/(c^2+4). '
   'Thus positive sign persistence is weaker than productive-cycle persistence: the cycle can remain positive while becoming exponentially sterile.'
 ),
 'rows':rows,'asymptotic':asym,
},indent=2,allow_nan=False))
