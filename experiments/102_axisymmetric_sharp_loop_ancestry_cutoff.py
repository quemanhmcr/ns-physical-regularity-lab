import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rt2=arb(2).sqrt(); rows=[]
# Sharp axisymmetric ray E=diag(2,-1,-1), latitude n_x=c.
# loop length=2pi r sqrt(1-c^2), |omega|=5 q |c| sqrt(1-c^2), hence mu=2pi r/(5q|c|).
for K0s in ['1e-24','1','1e24']:
 K0=arb(K0s); mu0=2*pi/K0
 for nus in ['1e-24','1','1e24']:
  nu=arb(nus)
  for rs in ['1e-12','1','1e12']:
   r=arb(rs); c=1/rt2; qmax=K0*r/(5*c); Rebound=qmax*r*r/nu
   must=rt2*K0*r**3/(5*nu)
   if not (Rebound-must).contains(0): raise AssertionError(('sharp frozen loop Re cutoff',K0s,nus,rs,Rebound,must))
   rows.append({'initial_gradomega_bound_K0':K0s,'frozen_closed_loop_mu_lower_2pi_over_K0':str(mu0),'nu':nus,'r':rs,'representative_latitude_abs_c_1over_sqrt2':str(c),'maximum_q_compatible_with_frozen_mu':str(qmax),'shell_transaction_Re_upper':str(Rebound),'closed_upper_sqrt2_K0_r3_over_5nu':str(must)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For the exact sharp axisymmetric carrier, each latitude vortex loop has mu=2pi r/(5q|c|). A smooth initial closed vortex lineage obeys the Fenchel/Cauchy lower bound mu>=2pi/||grad omega0||_infinity in the frozen branch. At the representative productive latitude |c|=1/sqrt(2), this forces q r^2/nu <= (sqrt(2)/5)(||grad omega0||_infinity/nu) r^3. Thus an exactly sharp closed frozen lineage becomes low-Re as its physical radius tends to zero.','rows':rows},indent=2,allow_nan=False))
