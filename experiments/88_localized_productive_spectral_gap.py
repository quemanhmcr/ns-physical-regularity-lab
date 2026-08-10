import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def certify_one(x,label,tol='1e-30'):
 t=arb(tol)
 if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
rows=[]
for p_int in [1,2,4,8,32]:
 p=arb(p_int)
 D=arb(1)/7-arb(2)/(p+7)+arb(1)/(2*p+7)
 N=p*p/(2*p+5)
 Lambda=N/D
 if not (Lambda>14): raise AssertionError(('weighted productive spectral gap violated',p_int,Lambda))
 for qs in ['1e-24','1','1e24']:
  q0=arb(qs)
  for Ls in ['1e-12','1','1e12']:
   L=arb(Ls)
   norm=q0*q0*L**3*D
   dirichlet=q0*q0*L*N
   rayleigh=dirichlet/norm
   certify_one((rayleigh*L*L)/Lambda,('Rayleigh quotient closed form',p_int,qs,Ls))
   if not (dirichlet >= (arb(14)/(L*L))*norm): raise AssertionError(('14/L2 coercivity',p_int,qs,Ls,dirichlet,norm))
   rows.append({'taper_power_p':p_int,'q0':qs,'L':Ls,'transaction_L2_radial_norm':str(norm),'radial_defect_Dirichlet_form':str(dirichlet),'dimensionless_Rayleigh_L2':str(rayleigh*L*L),'universal_lower_14':True})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For every regular localized productive profile q=r^2 k with q(L)=0, integration by parts gives <q,-L2 q>_{r^2dr}=integral r^6|k_prime|^2. The weighted Cauchy/Poincare estimate integral r^6 k^2 <= (L^2/14) integral r^6 k_prime^2 yields the universal source-scale gap -L2 >=14/L^2 on the localized toroidal l=2 sector. Polynomial tapers q=q0 x^2(1-x^p) all lie strictly above the bound.','rows':rows},indent=2,allow_nan=False))
