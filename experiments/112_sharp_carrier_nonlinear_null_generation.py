import json, os
from flint import arb, ctx
from fractions import Fraction as F
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rt2=arb(2).sqrt()
zero=arb(0); one=arb(1)

def padd(A,B):
 C=dict(A)
 for k,v in B.items(): C[k]=C.get(k,zero)+v
 return {k:v for k,v in C.items() if not v.contains(0) or v!=0}
def pscale(c,A): return {k:c*v for k,v in A.items()}
def pmul(A,B):
 C={}
 for e,u in A.items():
  for f,v in B.items():
   k=(e[0]+f[0],e[1]+f[1],e[2]+f[2]); C[k]=C.get(k,zero)+u*v
 return C
def pder(A,j):
 C={}
 for e,v in A.items():
  if e[j]:
   ee=list(e); c=ee[j]; ee[j]-=1; C[tuple(ee)]=C.get(tuple(ee),zero)+c*v
 return C
def prad(A,r): return {e:v*r**sum(e) for e,v in A.items()}
def vadd(a,b): return tuple(padd(a[i],b[i]) for i in range(3))
def vscale(c,a): return tuple(pscale(c,a[i]) for i in range(3))
def cross(a,b): return (padd(pmul(a[1],b[2]),pscale(-1,pmul(a[2],b[1]))),padd(pmul(a[2],b[0]),pscale(-1,pmul(a[0],b[2]))),padd(pmul(a[0],b[1]),pscale(-1,pmul(a[1],b[0]))))
def curl(a): return (padd(pder(a[2],1),pscale(-1,pder(a[1],2))),padd(pder(a[0],2),pscale(-1,pder(a[2],0))),padd(pder(a[1],0),pscale(-1,pder(a[0],1))))
def directional(v,a):
 out=[]
 for i in range(3):
  q={}
  for j in range(3): q=padd(q,pmul(v[j],pder(a[i],j)))
  out.append(q)
 return tuple(out)
def vdotpoly(a,b):
 q={}
 for i in range(3): q=padd(q,pmul(a[i],b[i]))
 return q

def odddf(n):
 if n<=0:return 1
 out=1
 while n>0: out*=n; n-=2
 return out
def sphere_avg_monomial(e):
 a,b,c=e
 if a%2 or b%2 or c%2:return F(0)
 aa,bb,cc=a//2,b//2,c//2; N=aa+bb+cc
 return F(odddf(2*aa-1)*oddd(2*bb-1)*oddd(2*cc-1),oddd(2*N+1))
def oddd(n): return odddf(n)
def sphere_avg(P):
 out=zero
 for e,v in P.items():
  f=sphere_avg_monomial(e)
  if f: out += v*arb(f.numerator)/f.denominator
 return out
def avg_vnorm2(v): return sphere_avg(vdotpoly(v,v))
def avg_vdot(a,b): return sphere_avg(vdotpoly(a,b))

X=({(1,0,0):one},{(0,1,0):one},{(0,0,1):one})
r2=padd(padd(pmul(X[0],X[0]),pmul(X[1],X[1])),pmul(X[2],X[2]))
# symmetric capacity stationary-lock S for s=1, L=1
c=arb(3)/(2*rt2)
S=((one,zero,c),(zero,one,c),(c,c,-arb(2)))
Sx=[]
for i in range(3):
 q={}
 for j in range(3): q=padd(q,pscale(S[i][j],X[j]))
 Sx.append(q)
Sx=tuple(Sx)
qpoly=vdotpoly(X,Sx)
# u=(1-5r2/3)Sx+(2/3)q x
u=[]
for i in range(3):
 term=padd(Sx[i],pscale(-arb(5)/3,pmul(r2,Sx[i])))
 term=padd(term,pscale(arb(2)/3,pmul(qpoly,X[i])))
 u.append(term)
u=tuple(u)
omega=curl(u)
R=vadd(directional(omega,u),vscale(-1,directional(u,omega)))

# verify exact carrier curl formula structurally by comparing to -(14/3)x cross Sx
expected=vscale(-arb(14)/3,cross(X,Sx))
for i in range(3):
 diff=padd(omega[i],pscale(-1,expected[i]))
 if any(not v.contains(0) for v in diff.values()): raise AssertionError(('carrier curl mismatch',i,diff))

n=X
rows=[]
for rs in ['1e-6','0.1','0.5','0.8','1']:
 r=arb(rs); Rr=tuple(prad(x,r) for x in R); wr=tuple(prad(x,r) for x in omega)
 nxR=cross(n,Rr)
 Q=[[zero for _ in range(3)] for _ in range(3)]
 for i in range(3):
  for j in range(3):
   qij=padd(pmul(n[i],nxR[j]),pmul(nxR[i],n[j])); Q[i][j]=arb(3)/2*sphere_avg(qij)
 trace=Q[0][0]+Q[1][1]+Q[2][2]
 if not trace.contains(0): raise AssertionError(('Q_R trace',rs,trace))
 # minimum productive vorticity corresponding to transaction tensor Q_R
 Qn=[]
 for i in range(3):
  q={}
  for j in range(3): q=padd(q,pscale(Q[i][j],n[j]))
  Qn.append(q)
 prod=vscale(-arb(5)/3,cross(n,tuple(Qn)))
 null=vadd(Rr,vscale(-1,prod))
 R2=avg_vnorm2(Rr); P2=avg_vnorm2(prod); N2=avg_vnorm2(null); crossPN=avg_vdot(null,prod)
 Q2=sum(Q[i][j]*Q[i][j] for i in range(3) for j in range(3))
 if not crossPN.contains(0): raise AssertionError(('nonlinear residual sharp/null lost orthogonality',rs,crossPN))
 if not (P2-(arb(5)/9)*Q2).contains(0): raise AssertionError(('sharp nonlinear residual floor',rs,P2,Q2))
 if not (R2-(P2+N2)).contains(0): raise AssertionError(('nonlinear residual Pythagoras',rs,R2,P2,N2))
 if r>=arb('0.5') and not (N2>0): raise AssertionError(('expected nonlinear null generation away from center',rs,N2))
 W2=avg_vnorm2(wr)
 rows.append({'r_over_L':rs,'nonlinear_residual_mean_square':str(R2),'productive_projection_mean_square':str(P2),'transaction_null_residual_mean_square':str(N2),'null_fraction_of_residual':str(N2/R2),'productive_fraction_of_residual':str(P2/R2),'carrier_vorticity_mean_square':str(W2),'null_generation_rate_over_s_times_carrier':str((N2/W2).sqrt()) if W2>0 else None,'Q_R_Frobenius_squared':str(Q2),'sharp_floor_5over9_Q2':str((arb(5)/9)*Q2),'raw_null_prod_cross':str(crossPN)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'The exact sharp tangent carrier of the symmetric stationary-lock tensor is not an invariant manifold of its own Euler vorticity dynamics. Projecting the nonlinear residual R=(omega.grad)u-(u.grad)omega shell-by-shell onto the sharp Hodge transaction projector gives an exact orthogonal split R=R_prod+R_null. The transaction-null component is strictly nonzero away from the center and is generated at order s times the carrier vorticity per strain time. Thus an exactly/near-sharp critical winding carrier requires continual external cancellation or replenishment of dynamically generated null modes; instantaneous projector efficiency alone is not dynamically persistent.','rows':rows},indent=2,allow_nan=False))
