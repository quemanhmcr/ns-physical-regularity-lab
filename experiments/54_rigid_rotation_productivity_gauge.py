import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def norm(a): return dot(a,a).sqrt()
def mv(A,v): return tuple(sum(A[i][j]*v[j] for j in range(3)) for i in range(3))
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x > 1-t and x < 1+t): raise AssertionError((label,x))

def geom(p,R,q):
    rp=norm(p); rq=norm(q); rr=norm(R); xi=tuple(v/rp for v in p); xj=tuple(v/rq for v in q); n=tuple(v/rr for v in R)
    D=dot(p,cross(R,q)); T=dot(xi,cross(n,xj)); alpha=dot(xi,n); beta=dot(xj,n); P=(T*alpha)*(-T*beta)
    return D,T,alpha,beta,P

p=(arb(1),arb(2),arb(-3)); R=(arb(4),arb(-5),arb(7)); q=(arb(-2),arb(6),arb(1))
base=geom(p,R,q)
angles=['1e-24','1','1000','1e12']
rows=[]
for ts in angles:
    th=arb(ts); c=th.cos(); s=th.sin()
    Q=((c,-s,arb(0)),(s,c,arb(0)),(arb(0),arb(0),arb(1)))
    pr=mv(Q,p); Rr=mv(Q,R); qr=mv(Q,q); g=geom(pr,Rr,qr)
    for k,name in enumerate(['D','T','alpha','beta','P']): certify_one(g[k]/base[k],('rotation invariant '+name,ts))
    rows.append({'angle':ts,'D_ratio':str(g[0]/base[0]),'T_ratio':str(g[1]/base[1]),'alpha_ratio':str(g[2]/base[2]),'beta_ratio':str(g[3]/base[3]),'P_ratio':str(g[4]/base[4])})

# Infinitesimal common skew mode: L=p.R is exactly unchanged; all squared magnitudes are unchanged.
ws=['1e-24','1','1e24']
skew=[]
for ss in ws:
    w=arb(ss); W=((arb(0),-w,arb(0)),(w,arb(0),arb(0)),(arb(0),arb(0),arb(0)))
    pd=mv(W,p); Rd=mv(W,R); qd=mv(W,q)
    Ldot=dot(pd,R)+dot(p,Rd)
    pdnorm=2*dot(p,pd); qdnorm=2*dot(q,qd); Rdnorm=2*dot(R,Rd)
    Ddot=dot(pd,cross(R,q))+dot(p,cross(Rd,q))+dot(p,cross(R,qd))
    for val,label in [(Ldot,'Ldot'),(pdnorm,'p norm'),(qdnorm,'q norm'),(Rdnorm,'R norm'),(Ddot,'Ddot')]:
        if not val.contains(0): raise AssertionError(('common skew not gauge',ss,label,val))
    skew.append({'omega':ss,'Ldot':str(Ldot),'p_norm2_dot':str(pdnorm),'q_norm2_dot':str(qdnorm),'R_norm2_dot':str(Rdnorm),'Ddot':str(Ddot)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','rotation_cases':len(rows),'skew_cases':len(skew),
 'interpretation':(
  'Common rigid rotation is an exact gauge mode of productive pair geometry: D, T, both longitudinal projections and the full cycle product P are invariant under SO(3). '
  'Infinitesimally a common skew velocity gradient changes neither p.R, the three magnitudes, nor the oriented pair cell D. '
  'Reset action must therefore measure material shape deformation after common rigid motion is removed, not total angular motion.'
 ),'rotations':rows,'skew':skew
},indent=2,allow_nan=False))
