import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

# Exact spatially linear incompressible Euler calibration:
# A=[[-a/2,-W/2,0],[W/2,-a/2,0],[0,0,a]], W'=aW.
# A'+A^2 is symmetric, so a quadratic pressure closes Euler exactly.
# Lambda=exp(a t): W~Lambda, axial material length~Lambda,
# transverse radii~Lambda^(-1/2), circulation and material volume are fixed.
# For an isolated material annular collar with fixed b/r,
# E_h~Lambda and tau_nu=(b-r)^2/nu~1/Lambda.

strain_rates=['1e-24','1e-12','1','1e12','1e24']
nus=['1e-24','1e-12','1','1e12','1e24']
r0s=['1e-24','1e-12','1','1e12','1e24']
gammas=['1e-18','1','1e18']
lambdas=['1','1e3','1e12','1e30','1e60']
ratio_b_over_r=arb(2)
ell0=arb('1.25')
rows=[]

for a_s in strain_rates:
  a=arb(a_s)
  for nu_s in nus:
    nu=arb(nu_s)
    for r0_s in r0s:
      r0=arb(r0_s)
      b0=ratio_b_over_r*r0
      delta0=b0-r0
      for G_s in gammas:
        G=arb(G_s)
        W0=G/(pi*r0*r0)
        Eh0=G*G*ell0/(4*pi)*(b0/r0).log()
        tau0=delta0*delta0/nu
        chi0=a*delta0*delta0/nu
        for L_s in lambdas:
          L=arb(L_s)
          sqrtL=L.sqrt()
          W=W0*L
          Wdot=a*W
          r=r0/sqrtL
          b=b0/sqrtL
          ell=ell0*L
          delta=b-r
          circulation=pi*r*r*W
          volume=pi*r*r*ell
          volume0=pi*r0*r0*ell0
          Eh=G*G*ell/(4*pi)*(b/r).log()
          tau=delta*delta/nu
          chi=a*delta*delta/nu

          tests={
            'W_amplification':W/W0/L,
            'axial_stretch':ell/ell0/L,
            'radius2_collapse':(r*r)/(r0*r0)*L,
            'circulation_conservation':circulation/G,
            'volume_conservation':volume/volume0,
            'collar_ratio_conservation':(b/r)/ratio_b_over_r,
            'harmonic_energy_amplification':Eh/Eh0/L,
            'viscous_clock_collapse':tau/tau0*L,
            'ancestry_clock_collapse':chi/chi0*L,
            'energy_clock_product':Eh*tau/(Eh0*tau0),
            'vorticity_ODE':Wdot/(a*W),
          }
          for name,val in tests.items():
            if not val.contains(1):
              raise AssertionError((name,a_s,nu_s,r0_s,G_s,L_s,val))

          # Euler closure: off-diagonal antisymmetric part of A'+A^2.
          # (A'+A^2)_12 = -Wdot/2 + a*W/2 = 0,
          # (A'+A^2)_21 =  Wdot/2 - a*W/2 = 0.
          e12=-Wdot/2+a*W/2
          e21=Wdot/2-a*W/2
          if not (e12.contains(0) and e21.contains(0)):
            raise AssertionError(('Euler matrix antisymmetric residual',a_s,G_s,L_s,e12,e21))

          rows.append({
            'a':a_s,'nu':nu_s,'r0':r0_s,'Gamma':G_s,'Lambda':L_s,
            'circulation_ratio':str(circulation/G),
            'volume_ratio':str(volume/volume0),
            'harmonic_energy_over_initial':str(Eh/Eh0),
            'viscous_clock_over_initial':str(tau/tau0),
            'ancestry_clock_Chi':str(chi),
            'Chi0':str(chi0),
            'energy_times_clock_ratio':str(Eh*tau/(Eh0*tau0)),
            'Euler_offdiag_residual_12':str(e12),
          })

# Anisotropic transverse attack: if lambda1*lambda2=1/Lambda,
# min(lambda1,lambda2) <= 1/sqrt(Lambda).  Parameterize lambda1=L^{-p},
# lambda2=L^{-(1-p)}, p in [0,1], and certify the shortest direction.
anis=[]
for L_s in ['1e3','1e12','1e30','1e60']:
  L=arb(L_s); iso=1/L.sqrt()
  for p_num in [0,1,2,3,4,5,6,7,8]:
    p=arb(p_num)/8
    # Use exp(-p log L) instead of noninteger power operator.
    l1=(-p*L.log()).exp()
    l2=(-(1-p)*L.log()).exp()
    prod=l1*l2*L
    if not prod.contains(1):
      raise AssertionError(('anisotropic area conservation',L_s,p_num,prod))
    # Choose the shorter factor from the exact rational exponent, not by
    # comparing potentially overlapping Arb intervals at p=1/2.
    shorter=l2 if p_num<=4 else l1
    gap=abs(p-arb('0.5'))
    expected_ratio=(-gap*L.log()).exp()
    ratio=shorter/iso
    if not (ratio/expected_ratio).contains(1):
      raise AssertionError(('anisotropic shortest-width formula mismatch',L_s,p_num,ratio,expected_ratio))
    if p_num==4:
      if not ratio.contains(1):
        raise AssertionError(('isotropic p=1/2 ratio must enclose 1',L_s,ratio))
    elif not (expected_ratio < 1):
      raise AssertionError(('anisotropic shortest direction should be strictly thinner',L_s,p_num,expected_ratio))
    anis.append({'Lambda':L_s,'p_eighths':p_num,
                 'lambda1':str(l1),'lambda2':str(l2),
                 'shorter_over_isotropic':str(ratio),
                 'area_ratio_certificate':str(prod)})

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'cases':len(rows),
  'anisotropy_cases':len(anis),
  'interpretation':'Exact linear Euler stretching of one frozen circulation ancestry amplifies vorticity and material length by Lambda while shrinking transverse area by 1/Lambda. In an isolated axisymmetric material collar the harmonic circulation energy grows by Lambda and its viscous crossing time falls by 1/Lambda, with circulation and material volume unchanged. Anisotropic transverse deformation cannot keep both widths larger than the Lambda^(-1/2) isotropic scale.',
  'anisotropy_attack':anis,
  'rows':rows,
},indent=2,allow_nan=False))
