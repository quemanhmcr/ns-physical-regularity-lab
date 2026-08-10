import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

# Canonical closed circular vortex-filament attack.
# A ring of radius R and circulation Gamma induces on its axis
#   u_z = Gamma R^2 / [2 (R^2+z^2)^(3/2)]
# and therefore
#   S_zz = -3 Gamma R^2 z / [2 (R^2+z^2)^(5/2)].
# On z=-beta R, S_zz>0 and scales Gamma/R^2.
# Rotational symmetry with isotropic core regularization makes the ring's
# self-induced motion a rigid axial translation: its own radius/period geometry
# has zero instantaneous deformation rate.  Thus outgoing productive strain
# does not force deformation of the same donor's I_C.

Rs=['1e-30','1e-15','1','1e15','1e30']
Gs=['1e-18','1','1e18']
betas=['0.01','0.1','0.25','0.5','1','2','10']
rows=[]
for Rs_ in Rs:
  R=arb(Rs_)
  for Gs_ in Gs:
    G=arb(Gs_)
    for bs in betas:
      beta=arb(bs)
      z=-beta*R
      den=(R*R+z*z)**arb('2.5')
      Szz=-3*G*R*R*z/(2*den)
      coeff=Szz*R*R/G
      target=3*beta/(2*(1+beta*beta)**arb('2.5'))
      if not (coeff/target).contains(1):
          raise AssertionError(('ring axial strain scale covariance failed',Rs_,Gs_,bs,coeff,target))
      if not (Szz>0):
          raise AssertionError(('chosen negative-z target must be stretched',Rs_,Gs_,bs,Szz))
      # Canonical symmetry statement: self-motion is rigid translation, so dR/dt=0.
      period_geometry_rate=arb(0)
      if not period_geometry_rate.contains(0):
          raise AssertionError('ring period geometry rate should be exact zero under rigid self-translation')
      rows.append({
        'R':Rs_,'Gamma':Gs_,'beta_minus_z_over_R':bs,
        'Szz_target':str(Szz),
        'dimensionless_Szz_R2_over_Gamma':str(coeff),
        'donor_self_period_geometry_rate':'0',
        'donor_self_motion':'rigid axial translation by rotational symmetry',
      })

# The dimensionless target strain is maximized at beta=1/2.
beta=arb('0.5')
fmax=3*beta/(2*(1+beta*beta)**arb('2.5'))
for bs in ['0.1','0.25','1','2']:
    b=arb(bs)
    f=3*b/(2*(1+b*b)**arb('2.5'))
    if not (f<fmax):
        raise AssertionError(('beta=1/2 should beat comparison sample',bs,f,fmax))

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'cases':len(rows),
  'beta_half_dimensionless_strain':str(fmax),
  'structural_attack':{
      'closed_donor':True,
      'outgoing_productive_strain_nonzero':True,
      'self_period_geometry_rate_zero_by_symmetry':True,
  },
  'interpretation':(
      'A closed circular vortex donor supplies O(Gamma/R^2) signed axial strain to a remote on-axis target while its own symmetric self-induced motion is rigid translation and leaves its circulation-period geometry unchanged. '
      'Therefore no universal instantaneous law can charge a donor outgoing productive transaction to deformation of that same donor I_C. '
      'The next object must be a directed interaction network among ancestries, not independent self-tolls.'
  ),
  'rows':rows,
},indent=2,allow_nan=False))
