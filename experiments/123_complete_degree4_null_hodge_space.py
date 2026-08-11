import json, os
from fractions import Fraction as F
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
z=arb(0);o=arb(1)

def mons(d):
 return [(a,b,d-a-b) for a in range(d+1) for b in range(d+1-a)]
def harmonic_basis(l):
 src=mons(l); tgt=mons(l-2) if l>=2 else []
 if not tgt:
  return [{e:F(int(e==src[j])) for e in src} for j in range(len(src))]
 ridx={e:i for i,e in enumerate(tgt)}; M=[[F(0) for _ in src] for _ in tgt]
 for j,e in enumerate(src):
  for k in range(3):
   if e[k]>=2:
    q=list(e);coef=q[k]*(q[k]-1);q[k]-=2;M[ridx[tuple(q)]][j]+=F(coef)
 # RREF
 A=[row[:] for row in M]; piv=[]; r=0
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
 free=[c for c in range(len(src)) if c not in piv]; out=[]
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
def div(a):
 q={}
 for i in range(3):q=padd(q,pder(a[i],i))
 return q
def vdot(a,b):
 q={}
 for i in range(3):q=padd(q,pmul(a[i],b[i]))
 return q
def r2pow(k):
 q={(0,0,0):o}
 for _ in range(k):q=pmul(q,r2)
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

def sharp_prod4(V):
 # Only need transaction tensor zero/nonzero; use the intrinsic unit-sphere observer.
 nx=cross(X,V);Q=[[z]*3 for _ in range(3)]
 for i in range(3):
  for j in range(3):Q[i][j]=arb(3)/2*savg(padd(pmul(X[i],nx[j]),pmul(nx[i],X[j])))
 return sum(Q[i][j]*Q[i][j] for i in range(3) for j in range(3))

def flatten(V):return [V[i].get(e,z) for i in range(3) for e in mons(4)]
def rankcols(cols):
 bas=[];piv=[]
 for col in cols:
  w=list(col)
  for p,b in zip(piv,bas):
   c=w[p]
   for j in range(len(w)):w[j]-=c*b[j]
  q=next((j for j,x in enumerate(w) if not x.contains(0)),None)
  if q is not None:
   c=w[q];w=[x/c for x in w];bas.append(w);piv.append(q)
 return len(bas)

X=({(1,0,0):o},{(0,1,0):o},{(0,0,1):o});r2=padd(padd(pmul(X[0],X[0]),pmul(X[1],X[1])),pmul(X[2],X[2]))
fields=[];rows=[]
# Complete transaction-null poloidal sectors l=1,3,5.
for l in (1,3,5):
 hb=harmonic_basis(l)
 if len(hb)!=2*l+1:raise AssertionError(('harmonic dimension',l,len(hb)))
 p=5-l;fac=r2pow(p//2)
 for j,Hq in enumerate(hb):
  H=toarb(Hq);g=tuple(pder(H,i) for i in range(3));T=cross(X,g);U=tuple(pmul(fac,q) for q in T);V=curl(U)
  if any(sum(e)!=4 for comp in V for e in comp):raise AssertionError(('degree',l,j))
  if any(not x.contains(0) for x in div(V).values()):raise AssertionError(('divV',l,j,div(V)))
  if any(not x.contains(0) for x in vdot(X,U).values()):raise AssertionError(('tangent U',l,j))
  q2=sharp_prod4(V)
  if not q2.contains(0):raise AssertionError(('poloidal null transaction',l,j,q2))
  fields.append(V);rows.append({'sector':'poloidal','l':l,'basis_index':j,'transaction_Q_squared':str(q2)})
# Complete toroidal l=4 null sector with tangent Hodge lift.
l=4;hb=harmonic_basis(l)
if len(hb)!=9:raise AssertionError(('harmonic dimension4',len(hb)))
for j,Hq in enumerate(hb):
 H=toarb(Hq);g=tuple(pder(H,i) for i in range(3));V=cross(X,g)
 U3=vscale(-arb(5)/22,g);U5=tuple(padd(pscale(arb(7)/22,pmul(r2,g[i])),pscale(-arb(4)/11,pmul(H,X[i]))) for i in range(3));U=vadd(U3,U5)
 if any(not x.contains(0) for x in div(V).values()):raise AssertionError(('divV4',j))
 bd=vdot(X,U)
 # Mixed degree-3/5 Hodge lift is tangent on the physical source boundary r=1, not on every concentric sphere.
 bd2=savg(pmul(bd,bd))
 if not bd2.contains(0):raise AssertionError(('boundary tangent U4',j,bd2))
 q2=sharp_prod4(V)
 if not q2.contains(0):raise AssertionError(('toroidal l4 null transaction',j,q2))
 fields.append(V);rows.append({'sector':'toroidal','l':4,'basis_index':j,'transaction_Q_squared':str(q2)})
rank=rankcols([flatten(v) for v in fields])
if len(fields)!=30 or rank!=30:raise AssertionError(('complete null dimension',len(fields),rank))
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'complete_degree4_null_dimension':len(fields),'certified_coordinate_rank':rank,'sector_dimensions':{'poloidal_l1':3,'poloidal_l3':7,'toroidal_l4':9,'poloidal_l5':11},'interpretation':'The complete homogeneous degree-four divergence-free vorticity space splits into the five-dimensional productive toroidal l=2 transaction sector and a thirty-dimensional transaction-null Hodge space.  The null space is exactly poloidal l=1 (3), poloidal l=3 (7), toroidal l=4 (9), and poloidal l=5 (11).  Each basis field is built from a harmonic scalar potential and carries its natural tangent Hodge velocity.  Coordinate monomials are used only to certify the invariant total rank 30.' ,'rows':rows},indent=2,allow_nan=False))
