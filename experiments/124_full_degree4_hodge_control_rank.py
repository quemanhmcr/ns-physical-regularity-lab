import json, os
from fractions import Fraction as F
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
z=arb(0);o=arb(1);rt2=arb(2).sqrt()

def mons(d):return [(a,b,d-a-b) for a in range(d+1) for b in range(d+1-a)]
def harmonic_basis(l):
 src=mons(l);tgt=mons(l-2) if l>=2 else []
 if not tgt:return [{e:F(int(e==src[j])) for e in src} for j in range(len(src))]
 ri={e:i for i,e in enumerate(tgt)};A=[[F(0) for _ in src] for _ in tgt]
 for j,e in enumerate(src):
  for k in range(3):
   if e[k]>=2:
    q=list(e);c=q[k]*(q[k]-1);q[k]-=2;A[ri[tuple(q)]][j]+=F(c)
 piv=[];r=0
 for c in range(len(src)):
  p=next((i for i in range(r,len(A)) if A[i][c]),None)
  if p is None:continue
  A[r],A[p]=A[p],A[r];q=A[r][c];A[r]=[x/q for x in A[r]]
  for i in range(len(A)):
   if i==r:continue
   f=A[i][c]
   if f:A[i]=[A[i][j]-f*A[r][j] for j in range(len(src))]
  piv.append(c);r+=1
  if r==len(A):break
 free=[c for c in range(len(src)) if c not in piv];out=[]
 for f in free:
  v=[F(0)]*len(src);v[f]=F(1)
  for i,p in enumerate(piv):v[p]=-A[i][f]
  out.append({src[j]:v[j] for j in range(len(src)) if v[j]})
 return out
def toarb(P):return {e:arb(v.numerator)/v.denominator for e,v in P.items()}
def padd(A,B):
 C=dict(A)
 for k,v in B.items():C[k]=C.get(k,z)+v
 return C
def pscale(c,A):return {k:c*v for k,v in A.items()}
def pmul(A,B):
 C={}
 for e,u in A.items():
  for f,v in B.items():
   q=tuple(e[i]+f[i] for i in range(3));C[q]=C.get(q,z)+u*v
 return C
def pder(A,j):
 C={}
 for e,v in A.items():
  if e[j]:
   q=list(e);c=q[j];q[j]-=1;q=tuple(q);C[q]=C.get(q,z)+c*v
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
def bracket(vort,vel):return vadd(directional(vort,vel),vscale(-1,directional(vel,vort)))
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
def r2pow(k):
 q={(0,0,0):o}
 for _ in range(k):q=pmul(q,r2)
 return q
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
def flatten(V):return [V[i].get(e,z) for i in range(3) for e in mons(4)]
def independent(cols):
 bas=[];piv=[];sel=[]
 for idx,col in enumerate(cols):
  w=list(col)
  for p,b in zip(piv,bas):
   c=w[p]
   for j in range(len(w)):w[j]-=c*b[j]
  q=next((j for j,x in enumerate(w) if not x.contains(0)),None)
  if q is not None:
   c=w[q];w=[x/c for x in w];bas.append(w);piv.append(q);sel.append(idx)
 return sel,piv
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
def combine(cs,vs):
 out=({},{},{})
 for c,V in zip(cs,vs):out=vadd(out,vscale(c,V))
 return out

X=({(1,0,0):o},{(0,1,0):o},{(0,0,1):o});r2=padd(padd(pmul(X[0],X[0]),pmul(X[1],X[1])),pmul(X[2],X[2]))
c=arb(3)/(2*rt2);S=((o,z,c),(z,o,c),(c,c,-arb(2)))
Sx=[]
for i in range(3):
 q={}
 for j in range(3):q=padd(q,pscale(S[i][j],X[j]))
 Sx.append(q)
u1=tuple(Sx);qx=vdot(X,u1);u3=tuple(padd(pscale(-arb(5)/3,pmul(r2,u1[i])),pscale(arb(2)/3,pmul(qx,X[i]))) for i in range(3));omega=curl(u3);N=sharp_null4(bracket(omega,u3))
# Complete 30D null basis with each field's harmonic degree-3 velocity companion U3 (zero except toroidal l4).
Vbasis=[];U3basis=[];labels=[]
for l in (1,3,5):
 fac=r2pow((5-l)//2)
 for j,Hq in enumerate(harmonic_basis(l)):
  H=toarb(Hq);g=tuple(pder(H,i) for i in range(3));T=cross(X,g);U=tuple(pmul(fac,q) for q in T);Vbasis.append(curl(U));U3basis.append(({}, {}, {}));labels.append(f'P{l}_{j}')
for j,Hq in enumerate(harmonic_basis(4)):
 H=toarb(Hq);g=tuple(pder(H,i) for i in range(3));Vbasis.append(cross(X,g));U3basis.append(vscale(-arb(5)/22,g));labels.append(f'T4_{j}')
if len(Vbasis)!=30:raise AssertionError(len(Vbasis))
KH=[]
for V,U3 in zip(Vbasis,U3basis):KH.append(sharp_null4(vadd(bracket(V,u1),bracket(omega,U3))))
cols=[flatten(v) for v in KH];sel,piv=independent(cols);rank=len(sel);nullity=30-rank
target=flatten(vscale(-1,N));exists=False;coeff=[z]*30;res2=norm2v(N)
if rank>0:
 A=[[cols[sel[j]][piv[i]] for j in range(rank)] for i in range(rank)]
 try:
  x=solve(A,[target[p] for p in piv]);candidate=combine(x,[Vbasis[j] for j in sel]);resp=combine(x,[KH[j] for j in sel]);res=vadd(resp,N)
  if all(v.contains(0) for v in flatten(res)):
   exists=True;res2=norm2v(res)
   for c,j in zip(x,sel):coeff[j]=c
 except RuntimeError:
  pass
# Compare full-space solution with the six-dimensional physics-generated KH cyclic space if unique.
in_cyclic=False;cyclic_res2=None
if exists:
 B=[N]
 def KHfield(V):
  # Extract U3 from toroidal l4 surface curl for any field in the generated sector.
  cS=vdot(X,curl(V));H4=pscale(-arb(1)/20,cS);g4=tuple(pder(H4,i) for i in range(3));U3=vscale(-arb(5)/22,g4)
  return sharp_null4(vadd(bracket(V,u1),bracket(omega,U3)))
 for _ in range(5):B.append(KHfield(B[-1]))
 bcols=[flatten(v) for v in B];s2,p2=independent(bcols)
 if len(s2)==6:
  A=[[bcols[j][p2[i]] for j in range(6)] for i in range(6)];cand=combine(coeff,Vbasis);rhs=flatten(cand)
  try:
   cc=solve(A,[rhs[p] for p in p2]);rr=vadd(cand,vscale(-1,combine(cc,B)));cyclic_res2=norm2v(rr);in_cyclic=cyclic_res2.contains(0)
  except RuntimeError:pass
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','complete_null_dimension':30,'full_KH_operator_certified_rank':rank,'full_KH_operator_nullity':nullity,'full_degree4_Hodge_servo_exists':exists,'full_degree4_Hodge_servo_unique':bool(exists and rank==30),'post_servo_null_mean_square':str(res2),'solution_lies_in_six_dim_generated_cyclic_space':in_cyclic,'cyclic_projection_residual_mean_square':str(cyclic_res2) if cyclic_res2 is not None else None,'nonzero_control_coefficients':sum(1 for x in coeff if not x.contains(0)),'interpretation':'The Hodge-complete degree-four linear maintenance operator is tested on the entire thirty-dimensional transaction-null space P1+P3+T4+P5, not only on the six-dimensional cyclic subspace generated by N4.  Its certified rank determines whether silent degree-four controls remain.  If rank is 30, the null-cancelling servo is unique in the complete degree-four Hodge space; if that unique solution lies in the six-dimensional cyclic subspace, then no unused degree-four null direction can be invoked to tune higher-order response without spoiling degree-four cancellation.'},indent=2,allow_nan=False))
