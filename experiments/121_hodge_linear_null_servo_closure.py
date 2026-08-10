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
def plap(A):
 C={}
 for j in range(3):C=padd(C,pder(pder(A,j),j))
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
def norm2s(P):return savg(pmul(P,P))
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
def hodge_U3(V):
 # Only the toroidal l=4 vorticity needs a lower-degree harmonic companion.
 cS=vdot(X,curl(V));H4=pscale(-arb(1)/20,cS);g4=tuple(pder(H4,i) for i in range(3))
 return vscale(-arb(5)/22,g4)
def KH(V):
 # Complete degree-four linear null response under the local Hodge lift.
 affine=vadd(mv(S,V),vscale(-1,directional(Sx,V)))
 U3=hodge_U3(V);harmonic_cross=vadd(directional(omega,U3),vscale(-1,directional(U3,omega)))
 return sharp_null4(vadd(affine,harmonic_cross))

mons=[]
for a in range(5):
 for b in range(5-a):mons.append((a,b,4-a-b))
coords=[(i,e) for i in range(3) for e in mons]
def flatten(V):return [V[i].get(e,z) for i,e in coords]
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
Sx=tuple(Sx);qx=vdot(X,Sx);u3=tuple(padd(pscale(-arb(5)/3,pmul(r2,Sx[i])),pscale(arb(2)/3,pmul(qx,X[i]))) for i in range(3));omega=curl(u3)
R4=vadd(directional(omega,u3),vscale(-1,directional(u3,omega)));N=sharp_null4(R4)
# Generate until rank stabilizes for several steps; permit the physics to choose the dimension.
vecs=[N];cols=[flatten(N)];rows=[];last_rank=0;stable=0
for k in range(14):
 if k>0:
  vecs.append(KH(vecs[-1]));cols.append(flatten(vecs[-1]))
 ps=pivots(cols);rank=len(ps)
 stable=stable+1 if rank==last_rank else 0;last_rank=rank
 rows.append({'iterate_k':k,'KH_vector_mean_square':str(norm2v(vecs[-1])),'certified_span_rank':rank})
# Use the final certified rank r and first r cyclic vectors to build the exact closure and solve KH V=-N if possible.
r=last_rank
B=vecs[:r];Br=KH(B[-1]);Bcols=[flatten(v) for v in B];ps=pivots(Bcols)
if len(ps)!=r:raise AssertionError(('rank basis failure',r,len(ps)))
Arel=[[Bcols[j][ps[i]] for j in range(r)] for i in range(r)];br=flatten(Br);arel=solve(Arel,[br[p] for p in ps]);res=vadd(Br,vscale(-1,combine(arel,B)))
if any(not x.contains(0) for x in flatten(res)):raise AssertionError(('KH closure relation failure',flatten(res)))
C=[[z for _ in range(r)] for _ in range(r)]
for j in range(r-1):C[j+1][j]=o
for i in range(r):C[i][r-1]=arel[i]
try:
 cs=solve(C,[-o]+[z]*(r-1));exists=True
except RuntimeError:
 cs=[];exists=False
if exists:
 V=combine(cs,B);post=vadd(KH(V),N);post2=norm2v(post)
 if not post2.contains(0):raise AssertionError(('Hodge linear servo residual',post2))
else:post2=norm2v(N)
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'certified_Hodge_linear_response_rank':r,'rank_stable_tail_steps':stable,'closure_relation_coefficients':[str(x) for x in arel],'full_Hodge_linear_null_servo_exists':exists,'servo_coefficients_in_KH_Krylov_basis':[str(x) for x in cs],'post_servo_degree4_null_residual_mean_square':str(post2),'interpretation':'The physically complete degree-four maintenance operator includes both affine transport of the null vorticity and the action of the unavoidable harmonic degree-three companion in its tangent Hodge velocity lift on the base l=2 vorticity.  The KH Krylov space is generated without a prescribed modal dimension.  Solving KH V=-N4 determines whether the first nonlinear null field can still be cancelled when the velocity-vorticity Hodge relation is respected.  Only after this closure is established is it meaningful to inspect degree-six/eight nonlinear responses.' ,'rows':rows},indent=2,allow_nan=False))
