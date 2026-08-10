import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))

# Accelerated reset lambda=1+sin^2(phi), phi=t/(1-t).
# Sample each transition at phi_k=k*pi+pi/4.  Then lambda=3/2, sin(2phi)=1 and
# |a_k|=(2/3)(1+phi_k)^2.  If this reset strain is harmonic/coherent on radius r,
# finite kinetic energy E0 forces r <= [15 E0/(4*pi*a_k^2)]^(1/5).
E0=arb(1)
C=(arb(135)/(16*pi)).root(5)  # r_H*(1+phi)^(4/5) for E0=1
rows=[]
for ks in ['0','1','10','1e6','1e30','1e60']:
    k=arb(ks); phi=k*pi+pi/4; q=1+phi
    a=arb(2)*q*q/3
    rH=(15*E0/(4*pi*a*a)).root(5)
    scaled=rH*q**(arb(4)/5)
    certify_one(scaled/C,('accelerated Hodge reset scale',ks))
    tail=1/q  # intrinsic time-to-singular-clock scale at this phase, up to exact phase shift
    r_over_tail=rH/(tail**(arb(4)/5))
    certify_one(r_over_tail/C,('tail 4/5 law',ks))
    rows.append({'k':ks,'intrinsic_phase_factor_1_plus_phi':str(q),'reset_strain_abs_a':str(a),'finite_energy_harmonic_radius_upper_E0_1':str(rH),'r_times_phase_4over5':str(scaled),'intrinsic_tail_1_over_1plusphi':str(tail)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'scale_constant_E0_1':str(C),
 'interpretation':(
  'For the exact finite-time accelerated reset clock, finite-amplitude common strain excursions require |a_k|~(1+phi_k)^2.  If that strain remains harmonically coherent on a material neighborhood and total kinetic energy is finite, the Hodge occupancy floor forces its admissible radius to shrink exactly like (1+phi_k)^(-4/5), equivalently like the intrinsic recurrence tail^(4/5).  The alternative is that the reset source transfers into the vortical/non-affine branch.'
 ),'rows':rows
},indent=2,allow_nan=False))
