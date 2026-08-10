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
def vadd(a,b):return tuple(padd(a[i],b[i]) for i in range(3))
def vscale(c,a):return tuple(pscale(c,a[i]) for i in range(3))
def cross(a,b):return (padd(pmul(a[1],b[2]),pscale(-1,pmul(a[2],b[1]))),padd(pmul(a[2],b[0]),pscale(-1,pmul(a[0],b[2]))),padd(pmul(a[0],b[1]),pscale(-1,pmul(a[1],b[0]))))
def curl(a):return (padd(pder(a[2],1),pscale(-1,pder(a[1],2))),padd(pder(a[0],2),pscale(-1,pder(a[2],0))),padd(pder(a[1],0),pscale(-1,pder(a[0],1))))
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
def norm2v(v):return savg(vdot(v,v))

def sharp_null4(V):
 nx=cross(X,V);Q=[[z]*3 for _ in range(3)]
 for i in range(3):
  for j in range(3):Q[i][j]=arb(3)/2*savg(padd(pmul(X[i],nx[j]),pmul(nx[i],X[j])))
 Qn=[]
 for i in range(3):
  q={}
  for j in range(3):q=padd(q,pscale(Q[i][j],X[j]))
  Qn.append(q)
 prod=vscale(-arb(5)/3,tuple(pmul(r2,c) for c in cross(X,tuple(Qn))))
 return vadd(V,vscale(-1,prod))
def K(V):return sharp_null4(vadd(mv(S,V),vscale(-1,directional(Sx,V))))

mons=[]
for a in range(5):
 for b in range(5-a):mons.append((a,b,4-a-b))
coords=[(i,e) for i in range(3) for e in mons]
def flatten(V):return [V[i].get(e,z) for i,e in coords]
def combine(coeffs,vecs):
 out=({},{},{})
 for c,V in zip(coeffs,vecs):out=vadd(out,vscale(c,V))
 return out
def pivots_from_columns(cols):
 basis=[];pivs=[]
 for col in cols:
  w=list(col)
  for p,b in zip(pivs,basis):
   c=w[p]
   for j in range(len(w)):w[j]-=c*b[j]
  pivot=next((j for j,x in enumerate(w) if not x.contains(0)),None)
  if pivot is not None:
   q=w[pivot];w=[x/q for x in w];basis.append(w);pivs.append(pivot)
 return pivs
def solve_square(A,b):
 n=len(A);M=[list(A[i])+[b[i]] for i in range(n)]
 for c in range(n):
  p=next((r for r in range(c,n) if not M[r][c].contains(0)),None)
  if p is None:raise RuntimeError(('singular_interval_matrix',c))
  M[c],M[p]=M[p],M[c];q=M[c][c]
  for j in range(c,n+1):M[c][j]/=q
  for r in range(n):
   if r==c:continue
   f=M[r][c]
   for j in range(c,n+1):M[r][j]-=f*M[c][j]
 return [M[i][n] for i in range(n)]

def coeff_residual(V):return flatten(V)

X=({(1,0,0):o},{(0,1,0):o},{(0,0,1):o});r2=padd(padd(pmul(X[0],X[0]),pmul(X[1],X[1])),pmul(X[2],X[2]))
c=arb(3)/(2*rt2);S=((o,z,c),(z,o,c),(c,c,-arb(2)))
Sx=[]
for i in range(3):
 q={}
 for j in range(3):q=padd(q,pscale(S[i][j],X[j]))
 Sx.append(q)
Sx=tuple(Sx);qx=vdot(X,Sx);u3=tuple(padd(pscale(-arb(5)/3,pmul(r2,Sx[i])),pscale(arb(2)/3,pmul(qx,X[i]))) for i in range(3))
omega=curl(u3);R4=vadd(directional(omega,u3),vscale(-1,directional(u3,omega)));N=sharp_null4(R4)
B=[N]
for _ in range(5):B.append(K(B[-1]))
B6=K(B[-1]);cols=[flatten(V) for V in B];pivs=pivots_from_columns(cols)
if len(pivs)!=6:raise AssertionError(('expected certified rank six',len(pivs),pivs))
Arel=[[cols[j][pivs[i]] for j in range(6)] for i in range(6)];brel=flatten(B6);arel=solve_square(Arel,[brel[p] for p in pivs])
relres=vadd(B6,vscale(-1,combine(arel,B)))
if any(not x.contains(0) for x in flatten(relres)):raise AssertionError(('Krylov closure relation not exact',flatten(relres)))
# Companion representation of K on the cyclic basis.
C=[[z for _ in range(6)] for _ in range(6)]
for j in range(5):C[j+1][j]=o
for i in range(6):C[i][5]=arel[i]
target=[-o,z,z,z,z,z]
try:
 servo_coeff=solve_square(C,target);servo_exists=True
except RuntimeError:
 servo_coeff=[];servo_exists=False
if servo_exists:
 Vservo=combine(servo_coeff,B);cancel=vadd(K(Vservo),N);flat=flatten(cancel)
 if any(not x.contains(0) for x in flat):raise AssertionError(('full affine null servo did not cancel N4',flat))
 servo_norm=norm2v(Vservo);cancel_norm=norm2v(cancel)
else:
 servo_norm=z;cancel_norm=norm2v(N)
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','krylov_rank':6,'closure_relation_K6_coefficients':[str(x) for x in arel],'full_affine_null_servo_exists':servo_exists,'servo_coefficients_in_physics_generated_Krylov_basis':[str(x) for x in servo_coeff],'servo_vorticity_mean_square':str(servo_norm),'post_servo_affine_null_residual_mean_square':str(cancel_norm),'interpretation':'The six-dimensional affine null-response space generated by K=P_null L_S is treated as the complete favorable servo space.  The exact Krylov closure relation determines the finite-dimensional action of K without choosing a modal ansatz.  Solving K V_servo=-N4 tests whether the entire first nonlinear null generation can be cancelled at affine order when all six physically generated servo directions are allowed.  Existence means affine maintenance is fully possible and the next obstruction, if any, must come from nonlinear interaction of the servo velocity with the base carrier.'},indent=2,allow_nan=False))
