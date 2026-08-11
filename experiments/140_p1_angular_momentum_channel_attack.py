import json, os
from fractions import Fraction as F
from flint import arb,ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160:raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import coupled46_hodge_core as H
import degree6_hodge_servo_core as C

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
 q=C.z
 for e,v in P.items():
  f=savgmono(e)
  if f:q+=v*arb(f.numerator)/f.denominator
 return q
def meanvec(V):return [savg(V[i]) for i in range(3)]
def meansq(V):
 m=meanvec(V);return sum((x*x for x in m),C.z),m

def r2pow(k,r2):
 q={(0,0,0):C.o}
 for _ in range(k):q=C.pmul(q,r2)
 return q
st=H.prepare()
# Nonsymmetric feedback states: no root or swap symmetry is imposed.
tests=[('zero',[0]*9),('basis_mix',[1,-2,3,4,-5,6,7,-8,9]),('asymmetric_large',[-7,1,11,-3,5,13,-2,17,-19])]
rows=[]
for name,vals in tests:
 y=[arb(v)/7 for v in vals];fb=H.feedback_map_native(st,y);hr=H.higher_responses_from_coupled(st,fb)
 for d in (8,10,12):
  N=hr[d][2];q,m=meansq(N)
  rows.append({'feedback_case':name,'degree':d,'sphere_mean_null_vorticity_vector':[str(v) for v in m],'sphere_mean_square':str(q),'P1_absent_by_mean_test':q.contains(0)})
# Sensitivity control: a genuine degree-eight P1 vorticity field has nonzero sphere mean.
a=({(0,0,0):C.o},{},{})
T=C.cross(st['X'],a);U=tuple(C.pmul(r2pow(4,st['r2']),q) for q in T);Vp1=C.curl(U);ctrl,ctrlm=meansq(Vp1)
if not (ctrl>0):raise AssertionError(('P1 sensitivity control',ctrl,ctrlm))
allzero=all(r['P1_absent_by_mean_test'] for r in rows)
if not allzero:raise AssertionError(('generic feedback generated P1',rows))
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'generic_nonsymmetric_feedback_P1_absent':allzero,'P1_sensitivity_control_mean_vector':[str(v) for v in ctrlm],'P1_sensitivity_control_mean_square':str(ctrl),'interpretation':'For a homogeneous divergence-free vorticity field, the poloidal l=1 sector is the only Hodge sector with nonzero sphere-mean vorticity; higher angular sectors have zero mean.  The observer is sensitive to an explicit degree-eight P1 control field.  Nevertheless the degree-eight, ten and twelve transaction-null emissions produced from several deliberately nonsymmetric degree-six feedback states all have zero sphere mean.  Thus the P1 absence seen at the certified coupled root is not caused merely by the x<->y symmetry or by root tuning.  It is a candidate structural selection rule of the tangent Hodge/Euler construction, naturally suggesting an angular-momentum/rotational conservation interpretation.'},indent=2,allow_nan=False))
