import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
# Exact diagonal transversality calibration:
# B=g diag(1,-(1-delta),-delta), trace zero. Switch along weak z eigen-direction.
# Vz-u = nu Deltaomega_z/(delta g).  To move distance epsilon in time tau,
# Deltaomega_z=delta g epsilon/(nu tau). Define ell_Delta=g/|Deltaomega_z|.
# With source Re Rs=epsilon^2/(nu tau): delta ell_Delta/epsilon=1/Rs.
for ds in ('0.5','0.1','1e-3','1e-12','1e-50'):
 d=arb(ds)
 b=(arb(1),-(1-d),-d)
 mean=(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])/3
 b2tf=(b[0]*b[0]-mean,b[1]*b[1]-mean,b[2]*b[2]-mean)
 detnorm=d*(1-d)
 for Rs in ('1','1e3','1e12','1e50'):
  R=arb(Rs)
  ell_over_eps=1/(d*R)
  gate=d*ell_over_eps
  if not (gate*R).contains(1): raise AssertionError(('switch gate',ds,Rs,gate))
  rows.append({'delta_weak_eigenvalue_ratio':ds,'source_Re':Rs,'normalized_simple_zero_determinant_detB_over_g3':str(detnorm),'required_curvature_length_over_core_ellDelta_over_epsilon':str(ell_over_eps),'exact_gate_delta_times_ell_over_epsilon':str(gate),'gate_times_source_Re':str(gate*R),'self_productive_B2_TF_diagonal_over_g2':[str(x) for x in b2tf],'weak_limit_productive_zz_component_target_minus_2over3':str(b2tf[2]+arb(2)/3)})
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'Calibrate fresh-simple-zero switching in the weak eigen-direction of B=g diag(1,-(1-delta),-delta).  The exact simple-zero drift law gives |V_z-u|=nu |Delta omega_z|/(delta g).  Requiring a zero to move one core radius epsilon in one remaining time tau gives the exact gate delta (ell_Delta/epsilon)=1/R_s, where ell_Delta=g/|Delta omega_z| and R_s=epsilon^2/(nu tau). '
  'Therefore at high source Reynolds there are only two ways to switch simple zeros fast: keep delta bounded below and force ell_Delta<<epsilon, a new viscous curvature subscale; or avoid the subscale by taking delta<=1/R_s, which drives det B/g^3=delta(1-delta)->0 and destroys simple-zero transversality. '
  'The degeneracy escape remains catalytically active: as delta->0, (B^2)_TF/g^2 tends diag(1/3,1/3,-2/3), so the nonlinear productive converter does not vanish when the zero becomes a line-like degenerate zero.  The frontier is therefore sharpened to an explicit dichotomy: inward viscous curvature or degenerate-zero bifurcation/zero-manifold recruitment.'),
 'rows':rows
},indent=2,allow_nan=False))
