import json, os, itertools
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def outer(a,b): return tuple(tuple(a[i]*b[j] for j in range(3)) for i in range(3))
def madd(*Ms): return tuple(tuple(sum(M[i][j] for M in Ms) for j in range(3)) for i in range(3))
def mscale(c,M): return tuple(tuple(c*M[i][j] for j in range(3)) for i in range(3))
def contract(A,B): return sum(A[i][j]*B[i][j] for i in range(3) for j in range(3))
def Aamp(u):
    return ((u[0]*u[0]-arb(1)/3,u[0]*u[1],u[0]*u[2]),(u[1]*u[0],u[1]*u[1]-arb(1)/3,u[1]*u[2]),(u[2]*u[0],u[2]*u[1],u[2]*u[2]-arb(1)/3))
def Mg(u,v,g): return madd(outer(u,v),outer(v,u),mscale(-g,outer(u,u)),mscale(-g,outer(v,v)))
def coord(M): return (M[0][0],M[1][1],M[0][1],M[0][2],M[1][2])
def parity(p):
    inv=sum(1 for i in range(len(p)) for j in range(i+1,len(p)) if p[i]>p[j])
    return -1 if inv%2 else 1
def det5(M):
    out=arb(0)
    for p in itertools.permutations(range(5)):
        term=arb(parity(p))
        for i in range(5): term*=M[i][p[i]]
        out+=term
    return out
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))

def make(alpha,beta,psi):
    sa=(1-alpha*alpha).sqrt(); sb=(1-beta*beta).sqrt()
    n=(arb(0),arb(0),arb(1)); a=(sa,arb(0),alpha)
    b=(sb*psi.cos(),-sb*psi.sin(),beta)
    return a,n,b

cases=[
 ('generic',None,None,None),
 ('capacity',1/arb(3).sqrt(),-1/arb(3).sqrt(),arb.pi()/2),
 ('near_coplanar',arb('0.4'),arb('-0.5'),arb('0.01')),
]
rows=[]
for name,alpha0,beta0,psi in cases:
    if name=='generic':
        a=(arb('0.8'),arb(0),arb('0.6')); n=(arb(0),arb(0),arb(1)); b=(arb('0.3'),arb('-0.4'),-arb(3).sqrt()/2)
    else:
        a,n,b=make(alpha0,beta0,psi)
    alpha=dot(a,n); beta=dot(b,n); gamma=dot(a,b); T=dot(a,cross(n,b))
    E=[Aamp(a),Aamp(b),Mg(a,n,alpha),Mg(b,n,beta),Mg(a,b,gamma)]
    C=[[coord(E[j])[i] for j in range(5)] for i in range(5)]
    Dc=det5(C); Dc_closed=T**4/3
    tol='1e-18' if name=='near_coplanar' else '1e-30'
    certify_one(Dc/Dc_closed,('pair strain frame coordinate volume',name),tol=tol)
    G=[[contract(E[i],E[j]) for j in range(5)] for i in range(5)]
    Dg=det5(G); Dg_closed=(arb(8)/3)*T**8
    certify_one(Dg/Dg_closed,('pair strain frame Frobenius Gram volume',name),tol=tol)
    if not (T*T>0): raise AssertionError(('noncoplanar case lost',name,T))
    rows.append({'case':name,'alpha':str(alpha),'beta':str(beta),'gamma':str(gamma),'T':str(T),'coordinate_frame_det':str(Dc),'T4_over_3':str(Dc_closed),'Frobenius_frame_Gram_det':str(Dg),'eight_over_3_T8':str(Dg_closed),'coordinate_ratio':str(Dc/Dc_closed),'Gram_ratio':str(Dg/Dg_closed)})

a,n,b=make(arb('0.4'),arb('-0.5'),arb(0)); T0=dot(a,cross(n,b))
if not T0.contains(0): raise AssertionError(('coplanar T',T0))
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'coplanar_structural_frame_volume':'0',
 'interpretation':(
  'A genuinely noncoplanar material pair generates five physical STF strain covectors: two magnitude-amplification tensors and three Gram-shape tensors.  Their invariant Frobenius Gram determinant is exactly (8/3) T^8, so they form a complete basis of incompressible symmetric strain space iff T is nonzero.  The pair therefore supplies its own physical strain frame and the frame loses rank exactly at coplanarity.'
 ),'rows':rows
},indent=2,allow_nan=False))
