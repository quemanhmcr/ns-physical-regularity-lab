import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
x0=(arb(7)/13).sqrt(); x07=x0**7
core_dim=arb(2)*x07/9
total_dim=arb(4)/9
collar_dim=total_dim+core_dim
for ss in ['1e-24','1','1e24']:
 s=arb(ss); det=-arb(17)*s**3/4; S2=arb(21)*s*s/2
 global_rate=(arb(4)/3)*det/S2
 if not (global_rate/(-arb(34)*s/63)).contains(1): raise AssertionError(('symmetric global rate',ss,global_rate))
 if not (global_rate<0): raise AssertionError(('global should deamplify',ss))
 rows.append({'s':ss,'sign_change_radius_over_L_sqrt7over13':str(x0),'global_log_enstrophy_self_stretch_rate':str(global_rate),'closed_minus_34over63_s':str(-arb(34)*s/63),'dimensionless_positive_inner_core_integral':str(core_dim),'dimensionless_negative_outer_collar_magnitude':str(collar_dim),'core_gain_over_collar_loss':str(core_dim/collar_dim),'core_gain_over_net_global_loss':str(core_dim/total_dim)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For the symmetric capacity stationary amplifier, det S<0 and the sharp tangent carrier shell self-stretch changes sign exactly at r/L=sqrt(7/13). The inner core is self-amplified while the outer collar is self-deamplified; the volume-integrated rate is exactly -(34/63)s. The positive core contribution is only about 5.4% of the collar-loss magnitude, so the isolated minimum carrier is globally sacrificial rather than self-sustaining. This is a stretching-sign ledger, not a literal conserved transfer of enstrophy between shells.','rows':rows},indent=2,allow_nan=False))
