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
zero_rows=[r for r in rows if r['feedback_case']=='zero']; nonsym=[r for r in rows if r['feedback_case']!='zero']
zero_absent=all(r['P1_absent_by_mean_test'] for r in zero_rows); nonsym_generated=all(not r['P1_absent_by_mean_test'] for r in nonsym)
if not zero_absent:raise AssertionError(('zero feedback unexpectedly generated P1',zero_rows))
if not nonsym_generated:raise AssertionError(('nonsymmetric P1 counterexample not separated from zero',nonsym))
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'zero_feedback_P1_absent':zero_absent,'nonsymmetric_feedback_generates_P1':nonsym_generated,'universal_P1_conservation_hypothesis_killed':True,'P1_sensitivity_control_mean_vector':[str(v) for v in ctrlm],'P1_sensitivity_control_mean_square':str(ctrl),'interpretation':'For a homogeneous divergence-free vorticity field, the poloidal l=1 sector is the only Hodge sector with nonzero sphere-mean vorticity; higher angular sectors have zero mean.  The observer is sensitive to an explicit degree-eight P1 control field.  The zero-feedback symmetric calibration has no P1 emission, but deliberately nonsymmetric degree-six feedback states generate large nonzero P1 sphere means at degrees eight, ten and twelve.  Therefore P1 is not universally protected by the tangent Hodge/Euler construction and the naive angular-momentum-conservation explanation is killed.  The P1 absence at the certified coupled root must instead come from the particular symmetry/compatibility of that equilibrium.'},indent=2,allow_nan=False))
