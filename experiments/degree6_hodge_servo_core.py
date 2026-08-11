from fractions import Fraction as F
from flint import arb
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
def bracket(vort,vel):return vadd(directional(vort,vel),vscale(-1,directional(vel,vort)))
def vdot(a,b):
 q={}
 for i in range(3):q=padd(q,pmul(a[i],b[i]))
 return q
def r2pow(k,r2):
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
def norm2v(v):return savg(vdot(v,v))
def norm2s(P):return savg(pmul(P,P))
def lapS(P,r2):return padd(pmul(r2,plap(P)),pscale(-1,padd(pD(pD(P)),pD(P))))
def spectral_project(P,l,allowed,r2):
 out=dict(P);lam=l*(l+1)
 for k in allowed:
  if k==l:continue
  lk=k*(k+1);out=pscale(arb(1)/(lk-lam),padd(lapS(out,r2),pscale(lk,out)))
 return out
def combine(cs,vs):
 out=({},{},{})
 for c,V in zip(cs,vs):out=vadd(out,vscale(c,V))
 return out
def flatten(V,d):return [V[i].get(e,z) for i in range(3) for e in mons(d)]
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

def setup_base():
 X=({(1,0,0):o},{(0,1,0):o},{(0,0,1):o});r2=padd(padd(pmul(X[0],X[0]),pmul(X[1],X[1])),pmul(X[2],X[2]))
 rt2=arb(2).sqrt();c=arb(3)/(2*rt2);S=((o,z,c),(z,o,c),(c,c,-arb(2)))
 Sx=[]
 for i in range(3):
  q={}
  for j in range(3):q=padd(q,pscale(S[i][j],X[j]))
  Sx.append(q)
 u1=tuple(Sx);qx=vdot(X,u1);u3=tuple(padd(pscale(-arb(5)/3,pmul(r2,u1[i])),pscale(arb(2)/3,pmul(qx,X[i]))) for i in range(3));omega=curl(u3)
 return X,r2,u1,u3,omega

def sharp_split(V,d,X,r2):
 nx=cross(X,V);Q=[[z]*3 for _ in range(3)]
 for i in range(3):
  for j in range(3):Q[i][j]=arb(3)/2*savg(padd(pmul(X[i],nx[j]),pmul(nx[i],X[j])))
 Qn=[]
 for i in range(3):
  q={}
  for j in range(3):q=padd(q,pscale(Q[i][j],X[j]))
  Qn.append(q)
 fac=r2pow((d-2)//2,r2);prod=vscale(-arb(5)/3,tuple(pmul(fac,c) for c in cross(X,tuple(Qn))))
 return prod,vadd(V,vscale(-1,prod))

def hodge_lift4(V,X,r2):
 h=vdot(X,V);H3=pscale(-arb(1)/216,plap(h));cS=vdot(X,curl(V));H4=pscale(-arb(1)/20,cS)
 g3=tuple(pder(H3,i) for i in range(3));UP5=tuple(pmul(r2,q) for q in cross(X,g3));g4=tuple(pder(H4,i) for i in range(3))
 U3=vscale(-arb(5)/22,g4);U5=tuple(padd(padd(pscale(arb(7)/22,pmul(r2,g4[i])),pscale(-arb(4)/11,pmul(H4,X[i]))),UP5[i]) for i in range(3))
 return U3,U5

def canonical_degree4_servo(X,r2,u1,u3,omega):
 R4base=bracket(omega,u3);N4=sharp_split(R4base,4,X,r2)[1]
 def KH(V):
  U3,_=hodge_lift4(V,X,r2);return sharp_split(vadd(bracket(V,u1),bracket(omega,U3)),4,X,r2)[1]
 B=[N4]
 for _ in range(5):B.append(KH(B[-1]))
 cols=[flatten(v,4) for v in B];sel,piv=independent(cols)
 if len(sel)!=6:raise AssertionError(('degree4 KH rank',len(sel)))
 Br=KH(B[-1]);A=[[cols[j][piv[i]] for j in range(6)] for i in range(6)];ar=solve(A,[flatten(Br,4)[p] for p in piv]);rel=vadd(Br,vscale(-1,combine(ar,B)))
 if any(not x.contains(0) for x in flatten(rel,4)):raise AssertionError('degree4 KH closure')
 C=[[z]*6 for _ in range(6)]
 for j in range(5):C[j+1][j]=o
 for i in range(6):C[i][5]=ar[i]
 cs=solve(C,[-o,z,z,z,z,z]);V4=combine(cs,B);U3,U5=hodge_lift4(V4,X,r2)
 R4=vadd(R4base,vadd(bracket(V4,u1),bracket(omega,U3)));N4post=sharp_split(R4,4,X,r2)[1]
 if not norm2v(N4post).contains(0):raise AssertionError(('degree4 servo residual',norm2v(N4post)))
 R6=vadd(vadd(bracket(V4,u3),bracket(omega,U5)),bracket(V4,U3));N6=sharp_split(R6,6,X,r2)[1]
 return V4,U3,U5,N6

def degree6_basis(X,r2):
 Vbasis=[];Ulow=[];Uhigh=[];labels=[];sectors=[]
 for l in (1,3,5,7):
  fac=r2pow((7-l)//2,r2)
  for j,Hq in enumerate(harmonic_basis(l)):
   H=toarb(Hq);g=tuple(pder(H,i) for i in range(3));T=cross(X,g);Uh=tuple(pmul(fac,q) for q in T);Vbasis.append(curl(Uh));Ulow.append(({}, {}, {}));Uhigh.append(Uh);labels.append(f'P{l}_{j}');sectors.append(('poloidal',l))
 for l in (4,6):
  qpow=6-l;facq=r2pow(qpow//2,r2);facqp2=r2pow((qpow+2)//2,r2);A=arb(qpow+l+3)/((qpow+2)*(qpow+2*l+3));B=-arb(l)/(qpow+2*l+3);C=-arb(l+1)/((qpow+2)*(qpow+2*l+3))
  for j,Hq in enumerate(harmonic_basis(l)):
   H=toarb(Hq);g=tuple(pder(H,i) for i in range(3));T=cross(X,g);V=tuple(pmul(facq,q) for q in T);Uh=tuple(padd(pscale(A,pmul(facqp2,g[i])),pscale(B,pmul(pmul(facq,H),X[i]))) for i in range(3));Ul=vscale(C,g)
   Vbasis.append(V);Ulow.append(Ul);Uhigh.append(Uh);labels.append(f'T{l}_{j}');sectors.append(('toroidal',l))
 return Vbasis,Ulow,Uhigh,labels,sectors

def degree6_diagonal_operator(V,Ul,sector,u1,omega,X,r2):
 R=bracket(V,u1)
 # Only the T6 harmonic degree-five companion contributes on the degree-six diagonal.
 if sector==('toroidal',6):R=vadd(R,bracket(omega,Ul))
 return sharp_split(R,6,X,r2)[1]

def degree6_lower_backreaction(Ul,sector,omega,X,r2):
 if sector==('toroidal',4):return sharp_split(bracket(omega,Ul),4,X,r2)[1]
 return ({},{},{})

def spectrum(V,d,X,r2):
 h=vdot(X,V);rad=tuple(pmul(h,X[i]) for i in range(3));T=vadd(V,vscale(-1,rad));divR=div(T);cor={}
 for i in range(3):
  for j in range(3):cor=padd(cor,pmul(X[i],pmul(X[j],pder(T[i],j))))
 dS=padd(divR,pscale(-1,cor));cS=vdot(X,curl(T));odd=list(range(1,d+2,2));even=list(range(2,d+1,2));pol={};tor={}
 for l in odd:
  q=spectral_project(dS,l,odd,r2);pol[l]=norm2s(q)/arb(l*(l+1))
 for l in even:
  q=spectral_project(cS,l,even,r2);tor[l]=norm2s(q)/arb(l*(l+1))
 return norm2s(h),pol,tor

def solve_degree6_servo():
 X,r2,u1,u3,omega=setup_base();V4,U3,U5,N6=canonical_degree4_servo(X,r2,u1,u3,omega);Vb,Ul,Uh,labels,sectors=degree6_basis(X,r2)
 K=[degree6_diagonal_operator(V,low,sec,u1,omega,X,r2) for V,low,sec in zip(Vb,Ul,sectors)];cols=[flatten(v,6) for v in K];sel,piv=independent(cols);rank=len(sel);coeff=[z]*len(Vb);exists=False;res=N6
 if rank:
  A=[[cols[sel[j]][piv[i]] for j in range(rank)] for i in range(rank)];target=flatten(vscale(-1,N6),6)
  try:
   x=solve(A,[target[p] for p in piv]);resp=combine(x,[K[j] for j in sel]);res=vadd(resp,N6)
   if all(v.contains(0) for v in flatten(res,6)):
    exists=True
    for c,j in zip(x,sel):coeff[j]=c
  except RuntimeError:pass
 V6=combine(coeff,Vb) if exists else ({},{},{})
 B4=combine(coeff,[degree6_lower_backreaction(low,sec,omega,X,r2) for low,sec in zip(Ul,sectors)]) if exists else ({},{},{})
 return {'X':X,'r2':r2,'u1':u1,'u3':u3,'omega':omega,'V4':V4,'N6':N6,'Vbasis':Vb,'Ulow':Ul,'Uhigh':Uh,'labels':labels,'sectors':sectors,'K66':K,'rank':rank,'coeff':coeff,'exists':exists,'residual':res,'V6':V6,'B4':B4}
