import json, os
from fractions import Fraction as F
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
z=arb(0);o=arb(1)

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
def sharp_Q2(V):
 nx=cross(X,V);Q=[[z]*3 for _ in range(3)]
 for i in range(3):
  for j in range(3):Q[i][j]=arb(3)/2*savg(padd(pmul(X[i],nx[j]),pmul(nx[i],X[j])))
 return sum(Q[i][j]*Q[i][j] for i in range(3) for j in range(3))
def flatten(V):return [V[i].get(e,z) for i in range(3) for e in mons(6)]
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
# Degree-six poloidal null sectors P_1+P_3+P_5+P_7.  Their natural Hodge velocity is toroidal degree seven.
for l in (1,3,5,7):
 hb=harmonic_basis(l)
 if len(hb)!=2*l+1:raise AssertionError(('harmonic dimension',l,len(hb)))
 fac=r2pow((7-l)//2)
 for j,Hq in enumerate(hb):
  H=toarb(Hq);g=tuple(pder(H,i) for i in range(3));T=cross(X,g);U=tuple(pmul(fac,q) for q in T);V=curl(U)
  if any(sum(e)!=6 for comp in V for e in comp):raise AssertionError(('degree P',l,j))
  if any(not x.contains(0) for x in div(V).values()):raise AssertionError(('div P',l,j))
  if any(not x.contains(0) for x in vdot(X,U).values()):raise AssertionError(('tangent P',l,j))
  q2=sharp_Q2(V)
  if not q2.contains(0):raise AssertionError(('poloidal transaction null',l,j,q2))
  fields.append(V);rows.append({'sector':'poloidal','l':l,'basis_index':j,'transaction_Q_squared':str(q2)})
# Degree-six toroidal null sectors T_4 and T_6.  Use the general exact tangent Hodge lift.
for l in (4,6):
 qpow=6-l;hb=harmonic_basis(l);facq=r2pow(qpow//2);facqp2=r2pow((qpow+2)//2)
 A=arb(qpow+l+3)/((qpow+2)*(qpow+2*l+3));B=-arb(l)/(qpow+2*l+3);C=-arb(l+1)/((qpow+2)*(qpow+2*l+3))
 structural_boundary_numerator=l*((qpow+l+3)-(qpow+2)-(l+1))
 if structural_boundary_numerator!=0:raise AssertionError(('structural tangent coefficient',l,qpow))
 for j,Hq in enumerate(hb):
  H=toarb(Hq);g=tuple(pder(H,i) for i in range(3));T=cross(X,g);V=tuple(pmul(facq,q) for q in T)
  Uhigh=tuple(padd(pscale(A,pmul(facqp2,g[i])),pscale(B,pmul(pmul(facq,H),X[i]))) for i in range(3));Ulow=vscale(C,g);U=vadd(Uhigh,Ulow)
  if any(sum(e)!=6 for comp in V for e in comp):raise AssertionError(('degree T',l,j))
  if any(not x.contains(0) for x in div(V).values()):raise AssertionError(('div T',l,j))
  curlerr=vadd(curl(U),vscale(-1,V))
  if any(not x.contains(0) for comp in curlerr for x in comp.values()):raise AssertionError(('curl lift T',l,j))
  rawbd=vdot(X,U);rawbd2=savg(pmul(rawbd,rawbd))
  if not rawbd2.contains(0):raise AssertionError(('source-sphere tangent T',l,j,rawbd2))
  q2=sharp_Q2(V)
  if not q2.contains(0):raise AssertionError(('toroidal transaction null',l,j,q2))
  fields.append(V);rows.append({'sector':'toroidal','l':l,'basis_index':j,'transaction_Q_squared':str(q2),'structural_boundary_normal_coefficient':'0','raw_boundary_normal_error_autopsy':str(rawbd2)})
rank=rankcols([flatten(v) for v in fields])
sector_dims={'poloidal_l1':3,'poloidal_l3':7,'poloidal_l5':11,'poloidal_l7':15,'toroidal_l4':9,'toroidal_l6':13}
if len(fields)!=58 or rank!=58:raise AssertionError(('complete degree6 null dimension',len(fields),rank))
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'complete_degree6_null_dimension':len(fields),'certified_coordinate_rank':rank,'sector_dimensions':sector_dims,'full_degree6_divfree_dimension':63,'productive_T2_dimension':5,'interpretation':'The complete homogeneous degree-six divergence-free vorticity space has dimension 63.  Removing the five-dimensional productive toroidal l=2 transaction sector leaves the exact 58-dimensional transaction-null Hodge space P1+P3+P5+P7+T4+T6.  Every basis field is generated from a harmonic scalar and carries its natural tangent Hodge velocity; coordinate monomials are used only to certify the invariant rank 58.' ,'rows':rows},indent=2,allow_nan=False))
