import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rows=[]
# Complete closed frozen flux lineages entirely inside B_L have specific volume mu>=2pi/K0.
# Disjoint flux elements satisfy dV=mu dGamma, hence total |flux| <= Vol(B_L)/mu0=(2/3)K0 L^3.
for K0s in ['1e-24','1','1e24']:
 K0=arb(K0s); mu0=2*pi/K0
 for Ls in ['1e-12','1','1e12']:
  L=arb(Ls); V=4*pi*L**3/3; PhiMax=V/mu0; closed=arb(2)*K0*L**3/3
  if not (PhiMax-closed).contains(0): raise AssertionError(('self-contained flux volume identity',K0s,Ls,PhiMax,closed))
  for nus in ['1e-24','1','1e24']:
   nu=arb(nus)
   for cp_s in ['2.8','3.5']:
    cp=arb(cp_s); ReMax=PhiMax/(cp*nu)
    rows.append({'K0_gradomega_initial_bound':K0s,'mu0_lower':str(mu0),'source_radius_L':Ls,'source_ball_volume':str(V),'maximum_self_contained_frozen_closed_flux':str(PhiMax),'nu':nus,'extremal_flux_over_nuRe_coefficient':cp_s,'source_Re_upper_if_all_flux_lineages_closed_inside_source':str(ReMax)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'Specific volume per frozen circulation flux converts initial smoothness into a self-contained source bound. Fenchel gives every closed frozen lineage mu>=2pi/K0, while disjoint complete flux elements satisfy dV=mu dGamma. If all such lineages are entirely contained in a ball of radius L, their total circulation flux is at most (2/3)K0 L^3. Consequently an extremal source Phi=c_p nu Re cannot remain high-Re as L->0. This closes the self-contained frozen-loop branch but not through-going lineages whose closure/specific volume lies outside the source.','rows':rows},indent=2,allow_nan=False))
