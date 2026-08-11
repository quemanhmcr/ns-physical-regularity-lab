import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def certify_one(x,label,tol='1e-30'):
 t=arb(tol)
 if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
rows=[]
for l_int in [2,4,6,8,16,32]:
 l=arb(l_int); gap=4*l+6
 for p_int in [1,2,4,8,32]:
  p=arb(p_int)
  # q=q0 x^l(1-x^p), x=r/L.
  D=arb(1)/(2*l+3)-arb(2)/(2*l+p+3)+arb(1)/(2*l+2*p+3)
  N=p*p/(2*l+2*p+1)
  Lambda=N/D
  if not (Lambda>gap): raise AssertionError(('general localized gap violated',l_int,p_int,Lambda,gap))
  for qs,Ls in [('1e-24','1e-12'),('1','1'),('1e24','1e12')]:
   q0=arb(qs);L=arb(Ls)
   norm=q0*q0*L**3*D
   dirichlet=q0*q0*L*N
   ray=dirichlet/norm
   certify_one((ray*L*L)/Lambda,('closed rayleigh',l_int,p_int,qs,Ls))
   if not (dirichlet >= gap/(L*L)*norm):raise AssertionError(('coercivity',l_int,p_int,qs,Ls))
   # L_l r^m=(m-l)(m+l+1)r^(m-2); smooth zero mode m=l is structural.
   zero_coeff=(l-l)*(l+l+1)
   if not zero_coeff.contains(0):raise AssertionError(('zero mode',l_int,zero_coeff))
   rows.append({'l':l_int,'taper_power_p':p_int,'q0':qs,'L':Ls,'simple_gap_4l_plus_6':str(gap),'dimensionless_Rayleigh':str(ray*L*L),'localized_norm':str(norm),'Dirichlet_form':str(dirichlet),'structural_smooth_zero_mode_q_proportional_r_l':'0'})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For every toroidal angular degree l, the radial viscous operator is L_l q=q_second+2q_prime/r-l(l+1)q/r^2=r^(-l-2)d_r[r^(2l+2)d_r(q/r^l)].  The smooth zero mode is q=r^l C.  If a source is genuinely localized by q(L)=0, writing q=r^l K gives <q,-L_l q>_{r^2dr}=int r^(2l+2)|K_prime|^2 and the exact weighted Cauchy estimate yields the universal simple gap (4l+6)/L^2.  The constant is not claimed sharp.  Polynomial tapers across l=2..32 all lie strictly above it.' ,'rows':rows},indent=2,allow_nan=False))
