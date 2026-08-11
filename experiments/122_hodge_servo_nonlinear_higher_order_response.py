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
def bracket(vort,vel):return vadd(directional(vort,vel),vscale(-1,directional(vel,vort)))
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
def lapS(P):return padd(pmul(r2,plap(P)),pscale(-1,padd(pD(pD(P)),pD(P))))
def spectral_project(P,l,allowed):
 out=dict(P);lam=l*(l+1)
 for k in allowed:
  if k==l:continue
  lk=k*(k+1);out=pscale(arb(1)/(lk-lam),padd(lapS(out),pscale(lk,out)))
 return out

def sharp_split(V,d):
 nx=cross(X,V);Q=[[z]*3 for _ in range(3)]
 for i in range(3):
  for j in range(3):Q[i][j]=arb(3)/2*savg(padd(pmul(X[i],nx[j]),pmul(nx[i],X[j])))
 Qn=[]
 for i in range(3):
  q={}
  for j in range(3):q=padd(q,pscale(Q[i][j],X[j]))
  Qn.append(q)
 fac={(0,0,0):o}
 for _ in range((d-2)//2):fac=pmul(fac,r2)
 prod=vscale(-arb(5)/3,tuple(pmul(fac,c) for c in cross(X,tuple(Qn))))
 return prod,vadd(V,vscale(-1,prod))
def hodge_lift4(V):
 h=vdot(X,V);H3=pscale(-arb(1)/216,plap(h));cS=vdot(X,curl(V));H4=pscale(-arb(1)/20,cS)
 g3=tuple(pder(H3,i) for i in range(3));UP5=tuple(pmul(r2,q) for q in cross(X,g3));g4=tuple(pder(H4,i) for i in range(3))
 U3=vscale(-arb(5)/22,g4);U5=[]
 for i in range(3):U5.append(padd(padd(pscale(arb(7)/22,pmul(r2,g4[i])),pscale(-arb(4)/11,pmul(H4,X[i]))),UP5[i]))
 return U3,tuple(U5)
def KH(V):
 U3,_=hodge_lift4(V);return sharp_split(vadd(bracket(V,u1),bracket(omega,U3)),4)[1]

def spectrum(V,d):
 h=vdot(X,V);rad=tuple(pmul(h,X[i]) for i in range(3));T=vadd(V,vscale(-1,rad));divR=div(T);cor={}
 for i in range(3):
  for j in range(3):cor=padd(cor,pmul(X[i],pmul(X[j],pder(T[i],j))))
 dS=padd(divR,pscale(-1,cor));cS=vdot(X,curl(T));odd=list(range(1,d+2,2));even=list(range(2,d+1,2));dp={l:spectral_project(dS,l,odd) for l in odd};cp={l:spectral_project(cS,l,even) for l in even}
 pol={l:norm2s(dp[l])/arb(l*(l+1)) for l in odd};tor={l:norm2s(cp[l])/arb(l*(l+1)) for l in even}
 radE=norm2s(h);tang=norm2v(T);polE=sum(pol.values(),z);torE=sum(tor.values(),z)
 if not (norm2v(V)-(radE+tang)).contains(0):raise AssertionError(('rad tan reconstruction',d))
 if not (tang-(polE+torE)).contains(0):raise AssertionError(('hodge reconstruction',d,tang,polE,torE))
 return radE,pol,tor

def mons_degree(d):
 out=[]
 for a in range(d+1):
  for b in range(d+1-a):out.append((a,b,d-a-b))
 return out
coords4=[(i,e) for i in range(3) for e in mons_degree(4)]
def flat4(V):return [V[i].get(e,z) for i,e in coords4]
def combine(cs,vs):
 out=({},{},{})
 for c,V in zip(cs,vs):out=vadd(out,vscale(c,V))
 return out
def pivots(cols):
 bas=[];ps=[]
 for col in cols:
  w=list(col)
  for p,b in zip(ps,bas):
   c=w[p]
   for j in range(len(w)):w[j]-=c*b[j]
  q=next((j for j,x in enumerate(w) if not x.contains(0)),None)
  if q is not None:
   c=w[q];w=[x/c for x in w];bas.append(w);ps.append(q)
 return ps
def solve(A,b):
 n=len(A);M=[list(A[i])+[b[i]] for i in range(n)]
 for c in range(n):
  p=next((r for r in range(c,n) if not M[r][c].contains(0)),None)
  if p is None:raise RuntimeError(('singular',c))
  M[c],M[p]=M[p],M[c];q=M[c][c]
  for j in range(c,n+1):M[c][j]/=q
  for r in range(n):
   if r==c:continue
   f=M[r][c]
   for j in range(c,n+1):M[r][j]-=f*M[c][j]
 return [M[i][n] for i in range(n)]

X=({(1,0,0):o},{(0,1,0):o},{(0,0,1):o});r2=padd(padd(pmul(X[0],X[0]),pmul(X[1],X[1])),pmul(X[2],X[2]))
c=arb(3)/(2*rt2);S=((o,z,c),(z,o,c),(c,c,-arb(2)))
Sx=[]
for i in range(3):
 q={}
 for j in range(3):q=padd(q,pscale(S[i][j],X[j]))
 Sx.append(q)
u1=tuple(Sx);qx=vdot(X,u1);u3=tuple(padd(pscale(-arb(5)/3,pmul(r2,u1[i])),pscale(arb(2)/3,pmul(qx,X[i]))) for i in range(3));omega=curl(u3)
R4base=bracket(omega,u3);N=sharp_split(R4base,4)[1]
# Build KH cyclic space until rank stabilizes at the validated six-dimensional value.
B=[N]
for _ in range(5):B.append(KH(B[-1]))
Br=KH(B[-1]);cols=[flat4(v) for v in B];ps=pivots(cols)
if len(ps)!=6:raise AssertionError(('KH rank not six',len(ps)))
Arel=[[cols[j][ps[i]] for j in range(6)] for i in range(6)];br=flat4(Br);arel=solve(Arel,[br[p] for p in ps]);rel=vadd(Br,vscale(-1,combine(arel,B)))
if any(not x.contains(0) for x in flat4(rel)):raise AssertionError(('KH relation fail',flat4(rel)))
C=[[z]*6 for _ in range(6)]
for j in range(5):C[j+1][j]=o
for i in range(6):C[i][5]=arel[i]
cs=solve(C,[-o,z,z,z,z,z]);V=combine(cs,B);U3,U5=hodge_lift4(V)
# Degree-four null cancellation check.
R4=vadd(R4base,vadd(bracket(V,u1),bracket(omega,U3)));_,N4post=sharp_split(R4,4)
if not norm2v(N4post).contains(0):raise AssertionError(('degree4 Hodge servo failed',norm2v(N4post)))
# Higher homogeneous responses of the actual Hodge servo.
R6=vadd(vadd(bracket(V,u3),bracket(omega,U5)),bracket(V,U3));P6,N6=sharp_split(R6,6)
R8=bracket(V,U5);P8,N8=sharp_split(R8,8)
rows=[]
for d,R,P,Nn in [(6,R6,P6,N6),(8,R8,P8,N8)]:
 rad,pol,tor=spectrum(Nn,d)
 row={'degree':d,'total_response_mean_square':str(norm2v(R)),'productive_projection_mean_square':str(norm2v(P)),'null_mean_square':str(norm2v(Nn)),'null_fraction':str(norm2v(Nn)/norm2v(R)) if norm2v(R)>0 else None,'radial_null_mean_square':str(rad)}
 for l,e in pol.items():row[f'poloidal_l{l}_energy']=str(e)
 for l,e in tor.items():row[f'toroidal_l{l}_energy']=str(e)
 rows.append(row)
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','degree4_post_servo_null_mean_square':str(norm2v(N4post)),'servo_vorticity_mean_square':str(norm2v(V)),'harmonic_U3_mean_square':str(norm2v(U3)),'vortical_U5_mean_square':str(norm2v(U5)),'interpretation':'Use the unique servo found inside the physics-generated six-dimensional KH space, lift it by the exact tangent Hodge velocity, and expand the complete Euler vorticity residual by physical homogeneity.  Degree four is certified null-free by construction.  The remaining degree-six response contains base-servo cross-interaction and servo interaction with its harmonic companion; degree eight is the servo self-interaction with its degree-five vortical velocity.  Sharp projection and intrinsic surface-Hodge spectra reveal which genuinely new response sectors are generated.  This is a canonical Hodge-closed servo calibration, not yet a universal no-control theorem because additional degree-four null controls outside the generated cyclic space have not been exhausted.','rows':rows},indent=2,allow_nan=False))
