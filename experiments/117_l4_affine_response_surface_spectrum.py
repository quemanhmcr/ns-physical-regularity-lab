import json, os
from flint import arb, ctx
from fractions import Fraction as F
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
z=arb(0);o=arb(1);rt2=arb(2).sqrt()

def padd(A,B):
 C=dict(A)
 for k,v in B.items():C[k]=C.get(k,z)+v
 return C
def pscale(c,A):return {k:c*v for k,v in A.items()}
def pmul(A,B):
 C={}
 for e,u in A.items():
  for f,v in B.items():
   k=tuple(e[i]+f[i] for i in range(3));C[k]=C.get(k,z)+u*v
 return C
def pder(A,j):
 C={}
 for e,v in A.items():
  if e[j]:
   ee=list(e);c=ee[j];ee[j]-=1;ee=tuple(ee);C[ee]=C.get(ee,z)+c*v
 return C
def pD(A):return {e:sum(e)*v for e,v in A.items()}
def plap(A):
 C={}
 for j in range(3):C=padd(C,pder(pder(A,j),j))
 return C
def vadd(a,b):return tuple(padd(a[i],b[i]) for i in range(3))
def vscale(c,a):return tuple(pscale(c,a[i]) for i in range(3))
def cross(a,b):return (padd(pmul(a[1],b[2]),pscale(-1,pmul(a[2],b[1]))),padd(pmul(a[2],b[0]),pscale(-1,pmul(a[0],b[2]))),padd(pmul(a[0],b[1]),pscale(-1,pmul(a[1],b[0]))))
def curl(a):return (padd(pder(a[2],1),pscale(-1,pder(a[1],2))),padd(pder(a[0],2),pscale(-1,pder(a[2],0))),padd(pder(a[1],0),pscale(-1,pder(a[0],1))))
def div(a):
 q={}
 for i in range(3):q=padd(q,pder(a[i],i))
 return q
def directional(v,a):
 out=[]
 for i in range(3):
  q={}
  for j in range(3):q=padd(q,pmul(v[j],pder(a[i],j)))
  out.append(q)
 return tuple(out)
def mv(S,v):
 out=[]
 for i in range(3):
  q={}
  for j in range(3):q=padd(q,pscale(S[i][j],v[j]))
  out.append(q)
 return tuple(out)
def vdot(a,b):
 q={}
 for i in range(3):q=padd(q,pmul(a[i],b[i]))
 return q
def odddf(n):
 if n<=0:return 1
 q=1
 while n>0:q*=n;n-=2
 return q
def savgmono(e):
 a,b,c=e
 if a%2 or b%2 or c%2:return F(0)
 aa,bb,cc=a//2,b//2,c//2;N=aa+bb+cc
 return F(odddf(2*aa-1)*odddf(2*bb-1)*odddf(2*cc-1),odddf(2*N+1))
def savg(P):
 q=z
 for e,v in P.items():
  f=savgmono(e)
  if f:q+=v*arb(f.numerator)/f.denominator
 return q
def norm2s(P):return savg(pmul(P,P))
def norm2v(v):return savg(vdot(v,v))
def inner(a,b):return savg(vdot(a,b))
def lapS(P):return padd(pmul(r2,plap(P)),pscale(-1,padd(pD(pD(P)),pD(P))))
def spectral_project(P,l,allowed):
 out=dict(P);lam=l*(l+1)
 for k in allowed:
  if k==l:continue
  lk=k*(k+1);out=pscale(arb(1)/(lk-lam),padd(lapS(out),pscale(lk,out)))
 return out

def sharp_split4(V):
 nx=cross(X,V);Q=[[z]*3 for _ in range(3)]
 for i in range(3):
  for j in range(3):Q[i][j]=arb(3)/2*savg(padd(pmul(X[i],nx[j]),pmul(nx[i],X[j])))
 Qn=[]
 for i in range(3):
  q={}
  for j in range(3):q=padd(q,pscale(Q[i][j],X[j]))
  Qn.append(q)
 prod=vscale(-arb(5)/3,tuple(pmul(r2,c) for c in cross(X,tuple(Qn))))
 return prod,vadd(V,vscale(-1,prod))
def spectrum(V):
 h=vdot(X,V);rad=tuple(pmul(h,X[i]) for i in range(3));T=vadd(V,vscale(-1,rad));divR=div(T);cor={}
 for i in range(3):
  for j in range(3):cor=padd(cor,pmul(X[i],pmul(X[j],pder(T[i],j))))
 dS=padd(divR,pscale(-1,cor));cS=vdot(X,curl(T));odd=[1,3,5];even=[2,4]
 dp={l:spectral_project(dS,l,odd) for l in odd};cp={l:spectral_project(cS,l,even) for l in even}
 pol=sum((norm2s(dp[l])/arb(l*(l+1)) for l in odd),z);tor=sum((norm2s(cp[l])/arb(l*(l+1)) for l in even),z)
 return {'total':norm2v(V),'radial':norm2s(h),'tangential':norm2v(T),'poloidal':pol,'toroidal':tor,'poloidal_l1':norm2s(dp[1])/2,'poloidal_l3':norm2s(dp[3])/12,'poloidal_l5':norm2s(dp[5])/30,'toroidal_l2':norm2s(cp[2])/6,'toroidal_l4':norm2s(cp[4])/20}

X=({(1,0,0):o},{(0,1,0):o},{(0,0,1):o});r2=padd(padd(pmul(X[0],X[0]),pmul(X[1],X[1])),pmul(X[2],X[2]))
c=arb(3)/(2*rt2);S=((o,z,c),(z,o,c),(c,c,-arb(2)))
Sx=[]
for i in range(3):
 q={}
 for j in range(3):q=padd(q,pscale(S[i][j],X[j]))
 Sx.append(q)
Sx=tuple(Sx);qx=vdot(X,Sx);u3=tuple(padd(pscale(-arb(5)/3,pmul(r2,Sx[i])),pscale(arb(2)/3,pmul(qx,X[i]))) for i in range(3))
omega=curl(u3);R4=vadd(directional(omega,u3),vscale(-1,directional(u3,omega)));_,N=sharp_split4(R4)
A=vadd(mv(S,N),vscale(-1,directional(Sx,N)));_,AN=sharp_split4(A);coef=inner(AN,N)/norm2v(N);orth=vadd(AN,vscale(-coef,N))
rows=[]
for name,V in [('generated_N4',N),('affine_null_response',AN),('orthogonal_new_null_direction',orth)]:
 sp=spectrum(V);rows.append({'field':name,**{k:str(v) for k,v in sp.items()}})
 # certify surface-Hodge reconstruction
 if not (sp['total']-(sp['radial']+sp['tangential'])).contains(0):raise AssertionError(('rad/tan',name,sp))
 if not (sp['tangential']-(sp['poloidal']+sp['toroidal'])).contains(0):raise AssertionError(('hodge',name,sp))
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'alignment_coefficient':str(coef),'interpretation':'The affine response of the generated l=4 null direction and the orthogonal servo remainder are resolved by the same intrinsic surface-Hodge operators.  This determines whether the one-channel failure stays inside the toroidal l=4 response sector or leaks immediately into radial/poloidal sectors.  If the latter remain structural zero, affine maintenance requires more angular directions but not a new angular degree; the true degree cascade must then come from nonlinear interaction with the localization field.','rows':rows},indent=2,allow_nan=False))
