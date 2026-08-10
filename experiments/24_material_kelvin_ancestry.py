import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

# Exact Lamb-Oseen material-circle calibration of viscous Kelvin circulation.
# Circle R=const is material because u is purely azimuthal.
# Gamma_R=Gamma(1-exp(-q)), q=R^2/(4 nu t).
# dGamma_R/dt = -Gamma*q*exp(-q)/t.
# curl(omega)_theta = omega_z*R/(2 nu t), so
# -nu integral_C curl(omega).dl gives the same quantity exactly.

Rs=['1e-30','1e-15','1','1e15','1e30']
nus=['1e-30','1','1e30']
Gs=['1e-20','1','1e20']
qs=['0.01','0.1','1','10','100']
rows=[]
for Rs_ in Rs:
  R=arb(Rs_)
  for nus_ in nus:
    nu=arb(nus_)
    for Gs_ in Gs:
      G=arb(Gs_)
      for qs_ in qs:
        q=arb(qs_)
        t=R*R/(4*nu*q)
        emq=(-q).exp()
        circ=G*(1-emq)
        direct=-G*q*emq/t
        omega=G/(4*pi*nu*t)*emq
        curlw_theta=omega*R/(2*nu*t)
        kelvin=-nu*(2*pi*R)*curlw_theta
        ratio=direct/kelvin
        if not ratio.contains(1):
            raise AssertionError(('Kelvin ancestry mismatch',Rs_,nus_,Gs_,qs_,direct,kelvin,ratio))
        # Dimensionless derivative removes arbitrary R,nu,G scales.
        reduced=(-direct)*t/G
        target=q*emq
        if not (reduced/target).contains(1):
            raise AssertionError(('reduced Kelvin defect lost scale covariance',Rs_,nus_,Gs_,qs_,reduced,target))
        if not (arb(0)<circ<G):
            raise AssertionError(('material circulation outside physical range',Rs_,nus_,Gs_,qs_,circ))
        rows.append({'R':Rs_,'nu':nus_,'Gamma':Gs_,'q':qs_,
                     'material_circulation':str(circ),
                     'direct_dGamma_dt':str(direct),
                     'viscous_kelvin_dGamma_dt':str(kelvin),
                     'identity_ratio':str(ratio),
                     'reduced_minus_dGamma_dt_t_over_Gamma':str(reduced)})

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'cases':len(rows),
  'interpretation':'For an exact Lamb-Oseen Navier-Stokes vortex, a fixed-radius circle is material and its explicit circulation derivative agrees identically with -nu times the line integral of curl(omega). Euler advection/stretching contributes no ancestry mutation. The identity remains scale-covariant over sixty decades of R and nu.',
  'rows':rows,
},indent=2,allow_nan=False))
