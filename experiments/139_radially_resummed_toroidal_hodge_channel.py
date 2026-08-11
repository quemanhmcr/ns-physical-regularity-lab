import json, os
from fractions import Fraction as F
from flint import arb,ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160:raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import degree6_hodge_servo_core as C
z=C.z;o=C.o

def re_z_power(l):
 # Re[(x+i y)^l], harmonic homogeneous scalar.
 P={}
 from math import comb
 for k in range(0,l+1,2):
  P[(l-k,k,0)]=arb(((-1)**(k//2))*comb(l,k))
 return P

def r2pow(k,r2):
 q={(0,0,0):o}
 for _ in range(k):q=C.pmul(q,r2)
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

def norm2v(V):return savg(C.vdot(V,V))
X=({(1,0,0):o},{(0,1,0):o},{(0,0,1):o});r2=C.padd(C.padd(C.pmul(X[0],X[0]),C.pmul(X[1],X[1])),C.pmul(X[2],X[2]))
profiles=[[(0,'1'),(2,'-0.3'),(4,'0.07')],[(0,'1e-20'),(2,'3e10'),(6,'-2e5'),(8,'0.125')],[(0,'1e20'),(4,'-2e-10'),(10,'3e-30')]]
rows=[]
for l in (2,4,6,8):
 H=re_z_power(l);g=tuple(C.pder(H,i) for i in range(3));T=C.cross(X,g)
 for pi,prof in enumerate(profiles):
  V=({},{},{});Uhi=({},{},{});C0=z;moment=z
  for q,cs in prof:
   c=arb(cs);rq=r2pow(q//2,r2);rq2=r2pow((q+2)//2,r2)
   V=C.vadd(V,tuple(C.pscale(c,C.pmul(rq,t)) for t in T))
   den=(q+2)*(q+2*l+3);A=c*arb(q+l+3)/den;B=-c*arb(l)/(q+2*l+3);cc=-c*arb(l+1)/den
   term=[]
   for i in range(3):term.append(C.padd(C.pscale(A,C.pmul(rq2,g[i])),C.pscale(B,C.pmul(rq,C.pmul(H,X[i])))))
   Uhi=C.vadd(Uhi,tuple(term));C0+=cc
   moment+=-arb(l+1)/arb(2*l+1)*c*(arb(1)/(q+2)-arb(1)/(q+2*l+3))
  U=C.vadd(Uhi,C.vscale(C0,g));curlerr=norm2v(C.vadd(C.curl(U),C.vscale(-1,V)));diverr=savg(C.pmul(C.div(U),C.div(U)));bd=C.vdot(X,U);bderr=savg(C.pmul(bd,bd))
  if not curlerr.contains(0):raise AssertionError(('curl',l,pi,curlerr))
  if not diverr.contains(0):raise AssertionError(('div',l,pi,diverr))
  if not bderr.contains(0):raise AssertionError(('boundary',l,pi,bderr))
  ratio=C0/moment if not moment.contains(0) else None
  if ratio is not None and not ratio.contains(1):raise AssertionError(('screened moment',l,pi,C0,moment,ratio))
  rows.append({'l':l,'profile_index':pi,'radial_profile_terms':[{'q':q,'coefficient':cs} for q,cs in prof],'harmonic_companion_C_l':str(C0),'screened_radial_moment':str(moment),'C_over_screened_moment':str(ratio) if ratio is not None else None,'curl_error':str(curlerr),'divergence_error':str(diverr),'boundary_tangency_error':str(bderr),'screen_kernel_exponent_2l_plus_1':2*l+1})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'A full toroidal angular Hodge channel omega=a(r) x cross grad H_l is radially resummed rather than treated as separate polynomial modes.  For a polynomial radial profile a(r), the exact tangent div-curl velocity equals a vortical radial part plus one harmonic companion C_l[a] grad H_l.  The companion is the screened radial moment C_l[a]=-(l+1)/(2l+1) integral_0^1 r a(r)[1-r^(2l+1)]dr.  Curl recovery, incompressibility, source-sphere tangency and the screened-moment identity are certified for mixed radial profiles and l=2,4,6,8.  At l=2 the screen is exactly 1-r^5, recovering the original Hodge transaction kernel.','rows':rows},indent=2,allow_nan=False))
