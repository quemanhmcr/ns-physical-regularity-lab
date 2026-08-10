import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def mv(A,v): return tuple(sum(A[i][j]*v[j] for j in range(3)) for i in range(3))
def add(a,b): return tuple(a[i]+b[i] for i in range(3))
def proj_perp(v,f):
    return tuple(f[i]-dot(v,f)*v[i] for i in range(3))
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x > 1-t and x < 1+t): raise AssertionError((label,x))

# Unit productive triad, not at the symmetric optimum.
a0=(arb('0.8'),arb(0),arb('0.6'))
n0=(arb(0),arb(0),arb(1))
b0=(arb('0.3'),arb('-0.4'),-arb(3).sqrt()/2)  # norm sqrt(.09+.16+.75)=1
for v,label in [(a0,'a'),(n0,'n'),(b0,'b')]:
    if not dot(v,v).contains(1): raise AssertionError(('unit vector',label,dot(v,v)))

# Residual shape drivers are represented directly in the common-spin quotient frame.
Sa=((arb('0.2'),arb('0.1'),arb('-0.05')),(arb('0.1'),arb('-0.15'),arb('0.07')),(arb('-0.05'),arb('0.07'),arb('-0.05')))
Sn=((arb('-0.1'),arb('0.04'),arb('0.03')),(arb('0.04'),arb('0.06'),arb('-0.02')),(arb('0.03'),arb('-0.02'),arb('0.04')))
Sb=((arb('0.05'),arb('-0.08'),arb('0.02')),(arb('-0.08'),arb('0.01'),arb('0.06')),(arb('0.02'),arb('0.06'),arb('-0.06')))
fa=proj_perp(a0,mv(Sa,a0)); fn=proj_perp(n0,mv(Sn,n0)); fb=proj_perp(b0,mv(Sb,b0))

alpha=dot(a0,n0); beta=dot(b0,n0); gamma=dot(a0,b0); T=dot(a0,cross(n0,b0))
F=1+2*alpha*beta*gamma-alpha*alpha-beta*beta-gamma*gamma; P=(-alpha*beta)*F
if not (P>0): raise AssertionError(('base pair not positive-product',P,alpha,beta,T))

ad=dot(fa,n0)+dot(a0,fn); bd=dot(fb,n0)+dot(b0,fn); gd=dot(fa,b0)+dot(a0,fb)
Fd=2*(ad*beta*gamma+alpha*bd*gamma+alpha*beta*gd-alpha*ad-beta*bd-gamma*gd)
Td=dot(fa,cross(n0,b0))+dot(a0,cross(fn,b0))+dot(a0,cross(n0,fb))
if not (Fd-2*T*Td).contains(0): raise AssertionError(('Gram determinant current vs triple product',Fd,2*T*Td))
Pd=-(ad*beta+alpha*bd)*F+(-alpha*beta)*Fd
Pd_direct=-2*T*Td*alpha*beta-T*T*(ad*beta+alpha*bd)
if not (Pd-Pd_direct).contains(0): raise AssertionError(('P current identity',Pd,Pd_direct))

# A common skew mode can be arbitrarily large; certify it separately as an exact gauge contribution.
rows=[]
for ws in ['1e-24','1','1e24']:
    w=arb(ws); W=((arb(0),-w,arb(0)),(w,arb(0),arb(0)),(arb(0),arb(0),arb(0)))
    Wa=mv(W,a0); Wn=mv(W,n0); Wb=mv(W,b0)
    ca=dot(Wa,n0)+dot(a0,Wn); cb=dot(Wb,n0)+dot(b0,Wn); cg=dot(Wa,b0)+dot(a0,Wb)
    cT=dot(Wa,cross(n0,b0))+dot(a0,cross(Wn,b0))+dot(a0,cross(n0,Wb))
    if not (ca.contains(0) and cb.contains(0) and cg.contains(0) and cT.contains(0)):
        raise AssertionError(('common spin leaked into Gram current',ws,ca,cb,cg,cT))
    rows.append({'common_spin_scale':ws,'alpha_common_rate':str(ca),'beta_common_rate':str(cb),'gamma_common_rate':str(cg),'T_common_rate':str(cT)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','common_spin_cases':len(rows),
 'base':{'alpha':str(alpha),'beta':str(beta),'gamma':str(gamma),'T':str(T),'Gram_det':str(F),'P':str(P)},
 'residual_current':{'alpha_dot':str(ad),'beta_dot':str(bd),'gamma_dot':str(gd),'Gram_det_dot':str(Fd),'two_T_Tdot':str(2*T*Td),'P_dot_gram':str(Pd),'P_dot_direct':str(Pd_direct)},
 'interpretation':(
  'The normalized productive state evolves through an intrinsic Gram current after common rigid spin is removed. '
  'The derivative of det G agrees exactly with 2 T Tdot from the residual triple-product current, and the resulting Pdot agrees with direct differentiation of -T^2 alpha beta. '
  'Common skew modes from 1e-24 to 1e24 contribute exactly zero to every Gram coordinate and to T.'
 ),'common_spin_gauge':rows
},indent=2,allow_nan=False))
