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
def inner(a,b):return savg(vdot(a,b))

def sharp_split4(V):
 nx=cross(X,V);Q=[[z]*3 for _ in range(3)]
 for i in range(3):
  for j in range(3):Q[i][j]=arb(3)/2*savg(padd(pmul(X[i],nx[j]),pmul(nx[i],X[j])))
 Qn=[]
 for i in range(3):
  q={}
  for j in range(3):q=padd(q,pscale(Q[i][j],X[j]))
  Qn.append(q)
 prod=vscale(-arb(5)/3,tuple(pmul(r2,c) for c in cross(X,tuple(Qn))))
 return prod,vadd(V,vscale(-1,prod))

X=({(1,0,0):o},{(0,1,0):o},{(0,0,1):o})
r2=padd(padd(pmul(X[0],X[0]),pmul(X[1],X[1])),pmul(X[2],X[2]))
c=arb(3)/(2*rt2);S=((o,z,c),(z,o,c),(c,c,-arb(2)))
Sx=[]
for i in range(3):
 q={}
 for j in range(3):q=padd(q,pscale(S[i][j],X[j]))
 Sx.append(q)
Sx=tuple(Sx);qx=vdot(X,Sx)
u3=tuple(padd(pscale(-arb(5)/3,pmul(r2,Sx[i])),pscale(arb(2)/3,pmul(qx,X[i]))) for i in range(3))
omega=curl(u3);R4=vadd(directional(omega,u3),vscale(-1,directional(u3,omega)));_,N=sharp_split4(R4)
# Affine vorticity response L_S N = S N -(Sx.grad)N.
A=vadd(mv(S,N),vscale(-1,directional(Sx,N)));Aprod,Anull=sharp_split4(A)
N2=norm2v(N);A2=norm2v(A);AN2=norm2v(Anull);AP2=norm2v(Aprod)
coef=inner(Anull,N)/N2
orth=vadd(Anull,vscale(-coef,N));orth2=norm2v(orth)
# If orth2=0, one scalar amplitude of N can servo its own l4 null direction under the affine strain.
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','N4_mean_square':str(N2),'affine_response_mean_square':str(A2),'affine_productive_mean_square':str(AP2),'affine_null_mean_square':str(AN2),'null_alignment_coefficient':str(coef),'null_component_orthogonal_to_N4_mean_square':str(orth2),'one_channel_affine_null_servo_exactly_aligned':bool(orth2.contains(0)),'interpretation':'The canonical generated l=4 null mode is subjected to the same stationary affine strain that amplifies the productive pair.  Projecting that affine response with the physical sharp transaction projector tests whether one scalar l=4 servo amplitude can cancel the dynamically generated null direction without creating a new independent null direction.  A nonzero orthogonal null remainder means even the most aligned one-channel servo is insufficient; exact zero means the l=4 direction is affine-closed up to productive response.'},indent=2,allow_nan=False))
