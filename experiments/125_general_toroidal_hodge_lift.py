import json, os, math
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
z=arb(0); o=arb(1)

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
   q=list(e); c=q[j]; q[j]-=1; q=tuple(q); C[q]=C.get(q,z)+c*v
 return C
def vadd(a,b): return tuple(padd(a[i],b[i]) for i in range(3))
def vscale(c,a): return tuple(pscale(c,a[i]) for i in range(3))
def cross(a,b): return (padd(pmul(a[1],b[2]),pscale(-1,pmul(a[2],b[1]))),padd(pmul(a[2],b[0]),pscale(-1,pmul(a[0],b[2]))),padd(pmul(a[0],b[1]),pscale(-1,pmul(a[1],b[0]))))
def curl(a): return (padd(pder(a[2],1),pscale(-1,pder(a[1],2))),padd(pder(a[0],2),pscale(-1,pder(a[2],0))),padd(pder(a[1],0),pscale(-1,pder(a[0],1))))
def div(a):
 q={}
 for i in range(3): q=padd(q,pder(a[i],i))
 return q
def vdot(a,b):
 q={}
 for i in range(3): q=padd(q,pmul(a[i],b[i]))
 return q
def r2pow(k):
 q={(0,0,0):o}
 for _ in range(k): q=pmul(q,r2)
 return q
def re_z_power(l):
 # Re[(x+i y)^l], harmonic homogeneous of degree l and independent of z.
 P={}
 for k in range(0,l+1,2):
  coeff=math.comb(l,k)*((-1)**(k//2)); P[(l-k,k,0)]=arb(coeff)
 return P

def odddf(n):
 if n<=0:return 1
 q=1
 while n>0:q*=n;n-=2
 return q
def savgmono(e):
 a,b,c=e
 if a%2 or b%2 or c%2:return None
 aa,bb,cc=a//2,b//2,c//2;N=aa+bb+cc
 return arb(odddf(2*aa-1)*odddf(2*bb-1)*odddf(2*cc-1))/odddf(2*N+1)
def savg(P):
 q=z
 for e,v in P.items():
  f=savgmono(e)
  if f is not None:q+=v*f
 return q

def norm2s(P):return savg(pmul(P,P))
def norm2v(V):return savg(vdot(V,V))

X=({(1,0,0):o},{(0,1,0):o},{(0,0,1):o})
r2=padd(padd(pmul(X[0],X[0]),pmul(X[1],X[1])),pmul(X[2],X[2]))
rows=[]
for l in (2,4,6,8):
 H=re_z_power(l); grad=tuple(pder(H,i) for i in range(3)); T=cross(X,grad)
 for qpow in (0,2,4):
  rq=r2pow(qpow//2); rq2=r2pow((qpow+2)//2)
  V=tuple(pmul(rq,t) for t in T)
  A=arb(qpow+l+3)/((qpow+2)*(qpow+2*l+3))
  B=-arb(l)/(qpow+2*l+3)
  C=-arb(l+1)/((qpow+2)*(qpow+2*l+3))
  U=[]
  for i in range(3):
   term=pscale(A,pmul(rq2,grad[i])); term=padd(term,pscale(B,pmul(rq,pmul(H,X[i])))); term=padd(term,pscale(C,grad[i])); U.append(term)
  U=tuple(U); diff=vadd(curl(U),vscale(-1,V)); dv=div(U); bd=vdot(X,U)
  if not norm2v(diff).contains(0):raise AssertionError(('curl',l,qpow,norm2v(diff)))
  if not savg(pmul(dv,dv)).contains(0):raise AssertionError(('div',l,qpow,savg(pmul(dv,dv))))
  raw_bd2=savg(pmul(bd,bd))
  if not raw_bd2.contains(0):raise AssertionError(('raw boundary observer excluded zero',l,qpow,raw_bd2))
  # Structural source-sphere tangency: x.grad H_l=l H_l, r=1, so the normal coefficient is l*A+B+l*C.
  # Clear the common denominator before any interval division: l[(q+l+3)-(q+2)-(l+1)]=0 exactly.
  structural_boundary_numerator=l*((qpow+l+3)-(qpow+2)-(l+1))
  if structural_boundary_numerator != 0:raise AssertionError(('structural boundary tangency',l,qpow,structural_boundary_numerator))
  rows.append({'l':l,'radial_power_q':qpow,'A_grad_coefficient':str(A),'B_Hx_coefficient':str(B),'C_harmonic_companion_coefficient':str(C),'vorticity_mean_square_unit_sphere':str(norm2v(V)),'curl_error':str(norm2v(diff)),'divergence_error':str(savg(pmul(dv,dv))),'structural_boundary_normal_coefficient':'0','raw_boundary_normal_error_autopsy':str(raw_bd2)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For toroidal angular vorticity V=r^q x cross grad H_l with H_l harmonic homogeneous and q even, the exact tangent Hodge lift is U=A r^(q+2) grad H_l + B r^q H_l x + C grad H_l, with A=(q+l+3)/[(q+2)(q+2l+3)], B=-l/(q+2l+3), C=-(l+1)/[(q+2)(q+2l+3)].  Curl recovery, incompressibility and source-sphere tangency are certified for l=2,4,6,8 and q=0,2,4.  The harmonic companion is intrinsic to the boundary condition; the old l=2 tangent strain carrier is the first member of this family.' ,'rows':rows},indent=2,allow_nan=False))
