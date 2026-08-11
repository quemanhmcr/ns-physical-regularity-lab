import json, os
from fractions import Fraction as F
from flint import arb,ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160:raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import coupled46_hodge_core as H
import degree6_hodge_servo_core as C
z=C.z;o=C.o;rt2=arb(2).sqrt()

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
def meanvec(V):return [savg(V[i]) for i in range(3)]
def ppow(P,n):
 q={(0,0,0):o}
 for _ in range(n):q=C.pmul(q,P)
 return q
def substitute(P,L):
 out={};cache=[{} for _ in range(3)]
 for i in range(3):cache[i][0]={(0,0,0):o}
 maxd=max((max(e) for e in P),default=0)
 for i in range(3):
  for k in range(1,maxd+1):cache[i][k]=C.pmul(cache[i][k-1],L[i])
 for e,v in P.items():
  q=C.pmul(C.pmul(cache[0][e[0]],cache[1][e[1]]),cache[2][e[2]])
  out=C.padd(out,C.pscale(v,q))
 return out
def transform_error(V,R,X):
 L=[]
 for i in range(3):
  q={}
  for j in range(3):q=C.padd(q,C.pscale(R[i][j],X[j]))
  L.append(q)
 VRx=[substitute(V[i],L) for i in range(3)];RV=[]
 for i in range(3):
  q={}
  for j in range(3):q=C.padd(q,C.pscale(R[i][j],V[j]))
  RV.append(q)
 return C.vadd(tuple(VRx),C.vscale(-1,tuple(RV)))
def rank3(A):
 cols=[[A[i][j] for i in range(len(A))] for j in range(3)]
 return len(C.independent(cols)[0])

st=H.prepare();sym=H.feedback_symmetry_basis(st);seq=C.solve_degree6_servo();y=[arb(seq['coeff'][i].mid()) for i in st['t4idx']];a,_=H.sym_coords(y,sym);a=[arb(v.mid()) for v in a]
for _ in range(8):
 g,_,_=H.reduced_feedback_native(st,sym,a);Jr=H.reduced_jacobian_native(st,sym,a);J=[[arb(Jr[i][j].mid()) for j in range(5)] for i in range(5)];dd=H.arbmat_solve(J,[-arb(v.mid()) for v in g]);a=[arb((a[i]+dd[i]).mid()) for i in range(5)]
_,fb,_=H.reduced_feedback_native(st,sym,a);hr=H.higher_responses_from_coupled(st,fb)
# R1: pi rotation about eigenvector (1,-1,0).  R2: pi rotation about the second eigenvector e+ +(sqrt2-1) ez.
R1=((z,-o,z),(-o,z,z),(z,z,-o))
t=rt2-o; den=o+t*t
# v=(1/sqrt2,1/sqrt2,t), projector vv^T/den.
v=(o/rt2,o/rt2,t)
R2=tuple(tuple(arb(2)*v[i]*v[j]/den-(o if i==j else z) for j in range(3)) for i in range(3))
# Both are proper involutive orthogonal symmetries of the stationary strain field.
fields={'u1':st['u1'],'u3':st['u3'],'omega2':st['omega'],'V4':fb['V4'],'V6':fb['V6'],'N8':hr[8][2],'N10':hr[10][2],'N12':hr[12][2]}
errs={}
for rn,R in [('R1',R1),('R2',R2)]:
 for name,V in fields.items():
  e=norm2v(transform_error(V,R,st['X']));errs[f'{rn}_{name}']=str(e)
  if not e.contains(0):raise AssertionError(('eigenframe symmetry',rn,name,e))
# The common fixed-vector space of R1 and R2 is zero, so any invariant sphere-mean vorticity must vanish.
A=[]
for R in (R1,R2):
 for i in range(3):A.append([R[i][j]-(o if i==j else z) for j in range(3)])
frank=rank3(A)
if frank!=3:raise AssertionError(('fixed-vector constraint rank',frank))
means={}
for d in (8,10,12):
 m=meanvec(hr[d][2]);q=sum((x*x for x in m),z);means[str(d)]={'mean':[str(x) for x in m],'mean_square':str(q)}
 if not q.contains(0):raise AssertionError(('P1 mean at symmetric root',d,q,m))
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','proper_eigenframe_rotations_tested':2,'common_fixed_vector_constraint_rank':frank,'common_fixed_axial_vector_dimension':0,'physical_field_symmetry_errors':errs,'higher_null_sphere_means':means,'interpretation':'The stationary strain has the proper D2 eigenframe symmetry of a symmetric tensor with distinct eigenvalues.  Two independent pi rotations are enough to eliminate every nonzero invariant vector.  The base carrier, the coupled V4/V6 servo, and the resulting N8/N10/N12 null emissions are directly certified invariant under both rotations.  Since a P1 vorticity component is detected by the sphere-mean vector, and no nonzero vector can be fixed by both rotations, P1 must vanish on this eigenframe-symmetric branch.  This explains the P1 absence without invoking a false universal angular-momentum conservation rule.'},indent=2,allow_nan=False))
