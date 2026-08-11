import json, os
from fractions import Fraction as F
from flint import arb,ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160:raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
z=arb(0)

def M(l,q):return -arb(l+1)/((q+2)*(q+2*l+3))
def Dcoef(l,q):
 # D_l r^q = q(q+2l+1) r^(q-2)
 return arb(q*(q+2*l+1))
profiles=[[(0,'1'),(2,'-0.3'),(4,'0.07')],[(0,'1e-20'),(2,'3e10'),(6,'-2e5'),(8,'0.125')],[(0,'1e20'),(4,'-2e-10'),(10,'3e-30')]]
rows=[]
for l in (2,4,6,8,12):
 for pi,prof in enumerate(profiles):
  CD=z;a0=z;a1=z
  for q,cs in prof:
   c=arb(cs);a1+=c
   if q==0:a0+=c
   if q>=2:CD+=c*Dcoef(l,q)*M(l,q-2)
  endpoint=arb(l+1)*(a0-a1);ratio=CD/endpoint if not endpoint.contains(0) else None
  if ratio is not None and not ratio.contains(1):raise AssertionError(('visc endpoint law',l,pi,CD,endpoint,ratio))
  rows.append({'l':l,'profile_index':pi,'C_l_of_D_l_a':str(CD),'endpoint_l_plus_1_a0_minus_aL':str(endpoint),'ratio':str(ratio) if ratio is not None else None})
 # Screened-null quadratic profile: silent statically, viscosity immediately makes feedback.
 cq=arb(2*(2*l+5))/arb(2*l+3);Cnull=M(l,0)-cq*M(l,2);visc=(l+1)*(arb(1)-(arb(1)-cq))
 if not Cnull.contains(0) or not (visc>0):raise AssertionError(('quadratic silent/viscous',l,Cnull,visc))
 rows.append({'l':l,'profile_index':'screened_null_quadratic','screened_moment':str(Cnull),'a0_minus_aL':str(cq),'C_l_of_D_l_a':str(visc),'viscosity_reactivates_feedback':True})
 # Doubly silent profile: C_l[a]=0 and a(0)=a(1), so first viscous derivative also zero.
 m0=M(l,0);m2=M(l,2);m4=M(l,4);A=-m0/(m2-m4);B=-A;Cd=m0+A*m2+B*m4;end=A+B
 first=arb(l+1)*(-end)
 # second feedback derivative C[D_l^2 a]=(l+1)(D_l a(0)-D_l a(1))=4(l+1)A(2l+5)
 second=arb(4*(l+1)*(2*l+5))*A
 if not Cd.contains(0) or not end.contains(0) or not first.contains(0) or second.contains(0):raise AssertionError(('double silent',l,Cd,end,first,second))
 rows.append({'l':l,'profile_index':'screened_and_first_viscous_silent','A_r2':str(A),'B_r4':str(B),'screened_moment':str(Cd),'a0_minus_aL':str(-end),'first_viscous_feedback_derivative_over_nu':str(first),'second_viscous_feedback_derivative_over_nu2':str(second),'viscous_silence_only_first_order':True})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For a toroidal angular Hodge channel omega=a(r) x cross grad H_l, the radial Laplacian is D_l a=a_double_prime+(2l+2)a_prime/r.  The screened Hodge feedback functional satisfies the exact adjoint identity C_l[D_l a]=(l+1)(a(0)-a(L)).  Thus pure viscosity changes the lower harmonic feedback moment only through center-to-source-boundary contrast.  A nonzero screened-null quadratic profile is immediately reactivated by viscosity.  One can additionally tune a quartic profile so C_l[a]=0 and a(0)=a(L), killing the first viscous feedback derivative, but its second derivative is nonzero; static lower-silence and even first-order viscous silence are therefore not the same as a viscosity-invariant hidden mode.','rows':rows},indent=2,allow_nan=False))
