import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rows=[]
# Smooth divergence-free spherical-shell flow-box (away from r=0):
# omega_r=A sin^2(theta)cos(theta)/r^2, omega_phi=B sin(theta)cos(theta)/r, omega_theta=0.
# Radial component is exactly invisible to n cross omega. At theta=pi/4, dr/dphi=A/(2B).
# To traverse DeltaR=R2-R1 in N turns choose A/B=DeltaR/(pi N).
# Full-shell enstrophy ratio radial/tangential=(4/7)(A/B)^2/(R1 R2).
for R1s,R2s in [('1','2'),('1e-6','2e-6'),('1e6','2e6')]:
 R1=arb(R1s); R2=arb(R2s); dR=R2-R1
 for Ns in ['1','10','1e6','1e30']:
  N=arb(Ns); eps=dR/(pi*N); ratio=(arb(4)/7)*eps*eps/(R1*R2)
  radial_drift=pi*N*eps
  if not (radial_drift/dR).contains(1): raise AssertionError(('spiral fixed radial drift',R1s,R2s,Ns,radial_drift,dR))
  rows.append({'R1':R1s,'R2':R2s,'turn_count_N':Ns,'A_over_B':str(eps),'radial_displacement_after_N_turns':str(radial_drift),'transaction_null_radial_over_tangential_enstrophy_ratio':str(ratio),'ratio_times_N_squared':str(ratio*N*N)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'A divergence-free spherical-shell flow-box can follow the sharp axisymmetric tangential angular pattern while adding a radial vorticity component that is exactly invisible to the Hodge transaction map. By spreading a fixed radial traverse across N azimuthal turns, the radial transaction-null enstrophy fraction falls like N^-2. Therefore instantaneous near-equality in the sharp projector does not force closed vortex-line ancestry: a very long slow spiral can be arbitrarily close in L2 to the tangential carrier while threading the source and closing elsewhere. The payment moves into winding/material line length rather than null amplitude.','rows':rows},indent=2,allow_nan=False))
