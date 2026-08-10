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
def pdegree(A,d): return {e:v for e,v in A.items() if sum(e)==d}
def vadd(a,b): return tuple(padd(a[i],b[i]) for i in range(3))
def vscale(c,a): return tuple(pscale(c,a[i]) for i in range(3))
def cross(a,b):
 return (padd(pmul(a[1],b[2]),pscale(-1,pmul(a[2],b[1]))),padd(pmul(a[2],b[0]),pscale(-1,pmul(a[0],b[2]))),padd(pmul(a[0],b[1]),pscale(-1,pmul(a[1],b[0]))))
def curl(a): return (padd(pder(a[2],1),pscale(-1,pder(a[1],2))),padd(pder(a[0],2),pscale(-1,pder(a[2],0))),padd(pder(a[1],0),pscale(-1,pder(a[0],1))))
def directional(v,a):
 out=[]
 for i in range(3):
  q={}
  for j in range(3): q=padd(q,pmul(v[j],pder(a[i],j)))
  out.append(q)
 return tuple(out)
def vdot(a,b):
 q={}
 for i in range(3): q=padd(q,pmul(a[i],b[i]))
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
def vnorm2(v): return savg(vdot(v,v))
def vinner(a,b): return savg(vdot(a,b))

def sharp_project(V,homdeg):
 # V is a homogeneous Cartesian vector polynomial of physical degree homdeg.
 # On radius r its transaction tensor scales r^homdeg.  At unit sphere n=x.
 n=X; nxV=cross(n,V)
 Q=[[z for _ in range(3)] for _ in range(3)]
 for i in range(3):
  for j in range(3):
   Q[i][j]=arb(3)/2*savg(padd(pmul(n[i],nxV[j]),pmul(nxV[i],n[j])))
 tr=sum(Q[i][i] for i in range(3))
 if not tr.contains(0): raise AssertionError(('transaction trace',homdeg,tr))
 Qn=[]
 for i in range(3):
  q={}
  for j in range(3): q=padd(q,pscale(Q[i][j],X[j]))
  Qn.append(q)
 # physical homogeneous extension of -(5/3)n x Q(r)n is -(5/3) r^(d-2) x x Qhat x.
 base=cross(X,tuple(Qn))
 power=homdeg-2
 if power<0 or power%2: raise AssertionError(('unexpected homogeneous degree',homdeg))
 fac={ (0,0,0):o }
 for _ in range(power//2): fac=pmul(fac,r2)
 prod=vscale(-arb(5)/3,tuple(pmul(fac,c) for c in base))
 null=vadd(V,vscale(-1,prod))
 V2=vnorm2(V); P2=vnorm2(prod); N2=vnorm2(null); PN=vinner(prod,null)
 Q2=sum(Q[i][j]*Q[i][j] for i in range(3) for j in range(3))
 if not PN.contains(0): raise AssertionError(('sharp orthogonality',homdeg,PN))
 if not (P2-(arb(5)/9)*Q2).contains(0): raise AssertionError(('sharp floor',homdeg,P2,Q2))
 if not (V2-(P2+N2)).contains(0): raise AssertionError(('pythagoras',homdeg,V2,P2,N2))
 return Q,prod,null,V2,P2,N2,PN

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
u1=Sx
u3=[]
for i in range(3):
 u3.append(padd(pscale(-arb(5)/3,pmul(r2,Sx[i])),pscale(arb(2)/3,pmul(qx,X[i]))))
u3=tuple(u3); omega=curl(u3)
R2=vadd(directional(omega,u1),vscale(-1,directional(u1,omega)))
R4=vadd(directional(omega,u3),vscale(-1,directional(u3,omega)))
# Structural polynomial-degree certificate.
for i in range(3):
 if any(sum(e)!=2 for e in R2[i]): raise AssertionError(('R2 degree leakage',i,R2[i]))
 if any(sum(e)!=4 for e in R4[i]): raise AssertionError(('R4 degree leakage',i,R4[i]))
Q2,P2,N2,V22,P22,N22,PN2=sharp_project(R2,2)
Q4,P4,N4,V24,P24,N24,PN4=sharp_project(R4,4)
# Degree-two affine response should be exactly on the productive transaction sector.
if not N22.contains(0): raise AssertionError(('affine degree-two response generated null sector',N22))
# Degree-four cubic self-interaction must generate a nonzero null sector.
if not (N24>0): raise AssertionError(('cubic degree-four response unexpectedly sharp',N24))
rows=[]
for d,V2v,P2v,N2v in [(2,V22,P22,N22),(4,V24,P24,N24)]:
 rows.append({'homogeneous_degree':d,'residual_mean_square_unit_sphere':str(V2v),'productive_projection_mean_square':str(P2v),'transaction_null_mean_square':str(N2v),'null_fraction':str(N2v/V2v) if V2v>0 else None})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For the exact sharp stationary-lock tangent carrier, the Euler vorticity residual splits structurally into a degree-two affine part and a degree-four cubic self-interaction.  The degree-two part lies exactly on the sharp productive transaction sector.  All first nonlinear departure from the sharp manifold is created by the homogeneous degree-four self-interaction, which has a strictly positive transaction-null component.  Thus null generation has a definite physical dilation order rather than being spread arbitrarily across the carrier.','rows':rows},indent=2,allow_nan=False))
