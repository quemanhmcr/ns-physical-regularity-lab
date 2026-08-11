import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); gamma=arb(4)/5; q=arb(2)/3; R0=arb(1); rows=[]
# Exact NS family u=(a(t)x+eps yz,-a(t)y,0), Delta u=0, arbitrary a(t).
# Set a=-d log lambda/dt and lambda=tau^(4/5)[1+sin^2 phi], phi=t/(1-t), tau=1/(1+phi).
# The material y-axis point (x,z)=(0,0), y0=q R0 follows y=lambda y0.
# Source radius R=R0 tau^(4/5), hence |y|/R=q[1+sin^2 phi], crossing 1 twice per pi-period.
for branch,phase0,sgn in [('outward',pi/4,1),('inward',3*pi/4,-1)]:
 for ks in ('0','1','10','1e6','1e30','1e60'):
  k=arb(ks); phi=k*pi+phase0; tau=1/(1+phi); f=arb(3)/2
  lam=tau**gamma*f; R=R0*tau**gamma; y=q*R0*lam
  ratio=abs(y)/R
  if not ratio.contains(1): raise AssertionError(('crossing ratio',branch,ks,ratio))
  # d/dt log(y/R)=sin(2phi) phi_dot/f = sgn*(2/3)/tau^2.
  dlog_ratio=arb(sgn)*arb(2)/(3*tau*tau)
  # common strain a=gamma/tau - sin(2phi)/(f tau^2).
  astrain=gamma/tau-arb(sgn)*arb(2)/(3*tau*tau)
  local_occupancy=astrain*astrain*R**5
  # dominant 4/5 law: a^2 R^5 tends 4/9 as tau->0 for R0=1.
  rows.append({
    'branch':branch,'k':ks,'phi':str(phi),'intrinsic_tail_tau':str(tau),'lambda':str(lam),'source_radius_R':str(R),
    'material_y':str(y),'abs_y_over_R':str(ratio),'crossing_log_speed':str(dlog_ratio),
    'common_strain_a':str(astrain),'local_Hodge_occupancy_scale_a2_R5':str(local_occupancy),
    'occupancy_over_asymptotic_4over9':str(local_occupancy/(arb(4)/9)),
  })
# Structural ancestry facts for this exact polynomial NS family:
# omega=(0,eps y,-eps z), curl omega=0, Delta omega=0.  Kelvin viscous current is identically zero.
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'source_shrink_exponent_gamma':'0.8','material_ratio_q_y0_over_R0':'0.6666666666666666666666666666666666666666666666667',
 'Kelvin_viscous_current_for_linear_vorticity':'0','material_flux_element_circulation_change_per_crossing':'0',
 'interpretation':(
  'Use the exact quadratic Navier-Stokes family u=(a x+eps yz,-a y,0), whose velocity and vorticity are harmonic polynomials so viscosity contributes zero locally.  Choose lambda=tau^(4/5)[1+sin^2 phi], phi=t/(1-t), a=-d log lambda/dt, and the shrinking source radius R=tau^(4/5). '
  'A material vorticity-line marker on the y-axis with y0=(2/3)R0 follows y=lambda y0, so |y|/R=(2/3)[1+sin^2 phi].  It crosses the shrinking source boundary exactly at phi=k pi+pi/4 and k pi+3pi/4 for every k, giving infinitely many alternating recrossings before t=1.  The same infinitesimal material circulation flux element is frozen because curl omega=0, so Kelvin mutation is exactly zero at every crossing. '
  'The accelerated crossing is locally compatible with the intrinsic harmonic energy scale: the common strain behaves like tau^-2 while R~tau^(4/5), hence a^2 R^5 tends the finite constant 4/9. '
  'Therefore infinite cumulative source-boundary circulation crossing does not imply infinite distinct ancestry or irreversible Kelvin renewal, even in an exact local Navier-Stokes family with the same 4/5 finite-energy scale.  This kills any attempt to use the module176 cumulative recruitment integral by itself as the contradiction. '
  'What remains stronger in module176 is the divergence of the instantaneous active cap flux.  Repeated sequential recrossing of one frozen lineage cannot by itself make the simultaneous flux inventory diverge; the next attack isolates multiplicity/packing at a fixed time.'),
 'rows':rows
},indent=2,allow_nan=False))
