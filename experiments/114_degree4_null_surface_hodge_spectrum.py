import json, os
from flint import arb, ctx
from fractions import Fraction as F
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
z=arb(0); o=arb(1); rt2=arb(2).sqrt()

def padd(A,B):
 C=dict(A)
 for k,v in B.items(): C[k]=C.get(k,z)+v
 return C
def pscale(c,A): return {k:c*v for k,v in A.items()}
def pmul(A,B):
 C={}
 for e,u in A.items():
  for f,v in B.items():
   k=tuple(e[i]+f[i] for i in range(3)); C[k]=C.get(k,z)+u*v
 return C
def pder(A,j):
 C={}
 for e,v in A.items():
  if e[j]:
   ee=list(e); c=ee[j]; ee[j]-=1; ee=tuple(ee); C[ee]=C.get(ee,z)+c*v
 return C
def pD(A): return {e:sum(e)*v for e,v in A.items()}
def plap(A):
 C={}
 for j in range(3): C=padd(C,pder(pder(A,j),j))
 return C
def vadd(a,b): return tuple(padd(a[i],b[i]) for i in range(3))
def vscale(c,a): return tuple(pscale(c,a[i]) for i in range(3))
def cross(a,b):
 return (padd(pmul(a[1],b[2]),pscale(-1,pmul(a[2],b[1]))),padd(pmul(a[2],b[0]),pscale(-1,pmul(a[0],b[2]))),padd(pmul(a[0],b[1]),pscale(-1,pmul(a[1],b[0]))))
def curl(a): return (padd(pder(a[2],1),pscale(-1,pder(a[1],2))),padd(pder(a[0],2),pscale(-1,pder(a[2],0))),padd(pder(a[1],0),pscale(-1,pder(a[0],1))))
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
 aa,bb,cc=a//2,b//2,c//2; N=aa+bb+cc
 return F(odddf(2*aa-1)*odddf(2*bb-1)*odddf(2*cc-1),odddf(2*N+1))
def savg(P):
 q=z
 for e,v in P.items():
  f=savgmono(e)
  if f:q+=v*arb(f.numerator)/f.denominator
 return q
def norm2s(P): return savg(pmul(P,P))
def norm2v(v): return savg(vdot(v,v))
def inner(a,b): return savg(vdot(a,b))
def lapS(P):
 # Delta_S = r^2 Delta_R3 - D^2 - D on the unit sphere.
 return padd(pmul(r2,plap(P)),pscale(-1,padd(pD(pD(P)),pD(P))))
def spectral_project(P,l,allowed):
 out=dict(P); lam=l*(l+1)
 for k in allowed:
  if k==l:continue
  lk=k*(k+1)
  out=pscale(arb(1)/(lk-lam),padd(lapS(out),pscale(lk,out)))
 return out

def sharp_project_degree4(V):
 n=X; nxV=cross(n,V)
 Q=[[z for _ in range(3)] for _ in range(3)]
 for i in range(3):
  for j in range(3): Q[i][j]=arb(3)/2*savg(padd(pmul(n[i],nxV[j]),pmul(nxV[i],n[j])))
 Qn=[]
 for i in range(3):
  q={}
  for j in range(3):q=padd(q,pscale(Q[i][j],X[j]))
  Qn.append(q)
 base=cross(X,tuple(Qn)); prod=vscale(-arb(5)/3,tuple(pmul(r2,c) for c in base))
 return prod,vadd(V,vscale(-1,prod))

X=({(1,0,0):o},{(0,1,0):o},{(0,0,1):o})
r2=padd(padd(pmul(X[0],X[0]),pmul(X[1],X[1])),pmul(X[2],X[2]))
c=arb(3)/(2*rt2)
S=((o,z,c),(z,o,c),(c,c,-arb(2)))
Sx=[]
for i in range(3):
 q={}
 for j in range(3):q=padd(q,pscale(S[i][j],X[j]))
 Sx.append(q)
Sx=tuple(Sx); qx=vdot(X,Sx)
u3=tuple(padd(pscale(-arb(5)/3,pmul(r2,Sx[i])),pscale(arb(2)/3,pmul(qx,X[i]))) for i in range(3))
omega=curl(u3)
R4=vadd(directional(omega,u3),vscale(-1,directional(u3,omega)))
prod,N=vsharp=sharp_project_degree4(R4)
# Surface-Hodge split of the degree-four null field on the unit source sphere.
h=vdot(X,N); radial=tuple(pmul(h,X[i]) for i in range(3)); T=vadd(N,vscale(-1,radial))
# Surface divergence is P_ij d_j T_i. For a tangent extension, div_S T = div T - n_i n_j d_j T_i.
divR=div(T)
cor={}
for i in range(3):
 for j in range(3): cor=padd(cor,pmul(X[i],pmul(X[j],pder(T[i],j))))
divS=padd(divR,pscale(-1,cor))
# Surface scalar curl is n dot curl T.
curlS=vdot(X,curl(T))
# The parity/degree structure permits radial and poloidal scalar sectors l=1,3,5; toroidal sector l=2,4.
odd=[1,3,5]; even=[2,4]
hparts={l:spectral_project(h,l,odd) for l in odd}
dparts={l:spectral_project(divS,l,odd) for l in odd}
cparts={l:spectral_project(curlS,l,even) for l in even}
# Reconstruct scalar observables on S2 and certify via L2 sphere norm.
def psum(parts):
 q={}
 for p in parts.values():q=padd(q,p)
 return q
for label,P,parts in [('radial',h,hparts),('divS',divS,dparts),('curlS',curlS,cparts)]:
 rem=padd(P,pscale(-1,psum(parts)))
 if not norm2s(rem).contains(0): raise AssertionError(('spectral reconstruction',label,norm2s(rem)))
# Productive l=2 was removed, so transaction-null tangential curl should have no l=2 sector.
if not norm2s(cparts[2]).contains(0): raise AssertionError(('l2 toroidal leakage survived sharp subtraction',norm2s(cparts[2])))
# Hodge energy reconstruction: tangential energy = sum ||div_l||^2/l(l+1)+||curl_l||^2/l(l+1).
radE=norm2s(h); tangE=norm2v(T); polE=z; torE=z
rows=[]
for l in odd:
 e=norm2s(dparts[l])/arb(l*(l+1)); polE+=e
 rows.append({'sector':'poloidal','l':l,'energy':str(e),'surface_divergence_L2':str(norm2s(dparts[l]))})
for l in even:
 e=norm2s(cparts[l])/arb(l*(l+1)); torE+=e
 rows.append({'sector':'toroidal','l':l,'energy':str(e),'surface_curl_L2':str(norm2s(cparts[l]))})
if not (tangE-(polE+torE)).contains(0): raise AssertionError(('surface Hodge energy reconstruction',tangE,polE,torE))
N2=norm2v(N)
if not (N2-(radE+tangE)).contains(0): raise AssertionError(('radial tangential split',N2,radE,tangE))
# Radial scalar spectrum is physically distinct and not counted in tangential Hodge energy.
for l in odd:
 rows.append({'sector':'radial','l':l,'energy':str(norm2s(hparts[l]))})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'degree4_null_total_mean_square':str(N2),'radial_mean_square':str(radE),'tangential_mean_square':str(tangE),'poloidal_mean_square':str(polE),'toroidal_mean_square':str(torE),'interpretation':'The first nonlinear transaction-null field generated by the sharp carrier has an intrinsic surface-Hodge spectrum.  After removing the productive toroidal l=2 transaction sector, the remaining degree-four null field splits into radial/poloidal odd spherical sectors and a higher toroidal sector.  The decomposition is obtained from the physical source sphere operators n dot, div_S, curl_S and Delta_S rather than from a chosen Cartesian modal basis.  This identifies which angular response channels a null-canceling network must actually control.','rows':rows},indent=2,allow_nan=False))
