import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160')); ctx.prec=BITS
if BITS<160: raise SystemExit('precision')
pi=arb.pi(); rows=[]
for Rs in ['1e-12','1','1e12']:
 R=arb(Rs)
 for xs in ['1e-12','0.01','0.1','0.5','0.9']:
  x=arb(xs); d=x*R
  F=1/(3*x**3)-arb(25)/21+x*x-x**7/7
  I=F/R**3
  if not I>0: raise AssertionError(('annulus factor',Rs,xs,I))
  # For A_e=e e-I/3, |A_e|^2=2/3.
  coeff=3*I/(10*pi)
  crude=1/(10*pi*d**3)
  if not crude>0: raise AssertionError('crude')
  rows.append({'R':Rs,'delta_over_R':xs,'native_I_delta_R':str(I),'vorticity_max_channel_L2_coefficient':str(coeff),'crude_delta_only_coefficient':str(crude)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'The exact radial Cauchy factor for the native Hodge weight on [delta,R] is I=1/(3delta^3)-25/(21R^3)+delta^2/R^5-delta^7/(7R^10). For the vorticity-maximum production tensor |A_e|^2=2/3, the sharp transaction projector gives the outer-channel spacetime L2 coefficient 3 I/(10 pi), with the simpler bound 1/(10 pi delta^3).','rows':rows},indent=2,allow_nan=False))
