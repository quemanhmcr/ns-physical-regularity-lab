import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

# Exact Burgers-vortex ancestry-clock calibration.
# Core scale delta_B^2=4 nu/a, hence a delta_B^2/nu=4 exactly.
# Compare this with a formal frozen-flux material collar delta_E^2=delta0^2/Lambda.
# The two scales meet at Lambda_cross=a delta0^2/(4 nu)=Chi0/4.
# This module does not claim the frozen-flux path is itself an NS solution;
# it verifies the scale-free crossover between exact Euler thinning kinematics
# and the exact steady NS stretch-diffuse balance.

strain_rates=['1e-30','1e-15','1','1e15','1e30']
nus=['1e-30','1e-15','1','1e15','1e30']
delta0s=['1e-30','1e-15','1','1e15','1e30']
gammas=['1e-18','1','1e18']
rows=[]

for a_s in strain_rates:
  a=arb(a_s)
  for nu_s in nus:
    nu=arb(nu_s)
    deltaB2=4*nu/a
    deltaB=deltaB2.sqrt()
    strain_time=1/a
    visc_clock_B=deltaB2/nu
    critical_clock=a*deltaB2/nu
    if not critical_clock.contains(4):
      raise AssertionError(('Burgers ancestry clock not pinned to 4',a_s,nu_s,critical_clock))
    if not (visc_clock_B/(4*strain_time)).contains(1):
      raise AssertionError(('Burgers viscous clock/strain time mismatch',a_s,nu_s,visc_clock_B,strain_time))

    for d0_s in delta0s:
      d0=arb(d0_s)
      chi0=a*d0*d0/nu
      Lcross=a*d0*d0/(4*nu)
      deltaE2=d0*d0/Lcross
      if not (deltaE2/deltaB2).contains(1):
        raise AssertionError(('Euler-thinning/Burgers crossover mismatch',a_s,nu_s,d0_s,deltaE2/deltaB2))
      if not (Lcross/(chi0/4)).contains(1):
        raise AssertionError(('crossover lost Chi0/4 scaling',a_s,nu_s,d0_s,Lcross,chi0))

      for G_s in gammas:
        G=arb(G_s)
        omega0=G*a/(4*pi*nu)
        # Exact Burgers central vorticity relative to strain is controlled only by Re_Gamma.
        reduced=omega0/a
        expected=G/(4*pi*nu)
        if not (reduced/expected).contains(1):
          raise AssertionError(('Burgers central-vorticity circulation-Re scaling failed',a_s,nu_s,G_s,reduced,expected))
        # Exact one-strain-time vorticity dissipation toll per axial length.
        toll=(G*G*a/(8*pi))*strain_time
        target=G*G/(8*pi)
        if not (toll/target).contains(1):
          raise AssertionError(('Burgers Gamma^2 toll mismatch',a_s,nu_s,G_s,toll,target))
        rows.append({
          'a':a_s,'nu':nu_s,'delta0':d0_s,'Gamma':G_s,
          'delta_B2':str(deltaB2),
          'Burgers_Chi':str(critical_clock),
          'viscous_clock_over_strain_time':str(visc_clock_B/strain_time),
          'Chi0':str(chi0),
          'Lambda_cross':str(Lcross),
          'deltaE2_over_deltaB2_at_cross':str(deltaE2/deltaB2),
          'central_omega_over_strain':str(reduced),
          'one_strain_toll_over_Gamma2_over_8pi':str(toll/target),
        })

# Explicit amplification clock: a frozen material collar with finite Chi0
# reaches Chi=1 after Lambda=Chi0.  At the Burgers Gaussian scale it reaches Chi=4.
clock=[]
for chi0_s in ['1e-12','1','4','1e3','1e12','1e30','1e60']:
  chi0=arb(chi0_s)
  for L_s in ['1','4','1e3','1e12','1e30','1e60']:
    L=arb(L_s)
    chi=chi0/L
    clock.append({'Chi0':chi0_s,'Lambda':L_s,'Chi_after_amplification':str(chi)})

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'cases':len(rows),
  'clock_cases':len(clock),
  'interpretation':'The exact steady Burgers vortex pins the transverse ancestry clock a*delta^2/nu to 4. Formal frozen-flux thinning from an initial collar delta0 reaches that exact NS stretch-diffuse scale at Lambda=Chi0/4, independent of absolute length scale. At the same balance the classic one-strain-time viscous vorticity toll is Gamma^2/(8*pi). This calibrates the persistence-to-diffusion crossover; it is not yet a universal 3D leakage theorem.',
  'amplification_clock':clock,
  'rows':rows,
},indent=2,allow_nan=False))
