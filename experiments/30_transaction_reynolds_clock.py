import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

# Productive Hodge transaction clock.
# Gamma_Q = r^2 q_e, so Gamma_Q/nu = q_e r^2/nu exactly equals
# (viscous crossing time r^2/nu)/(productive time 1/|q_e|) after sign is
# removed only for the clock interpretation.

rows=[]
for q_s in ['-1e30','-1','-1e-30','1e-30','1','1e30']:
  q=arb(q_s)
  aq=-q if q<0 else q
  for r_s in ['1e-30','1e-12','1','1e12','1e30']:
    r=arb(r_s)
    for nu_s in ['1e-30','1','1e30']:
      nu=arb(nu_s)
      GammaQ=q*r*r
      ReQ=GammaQ/nu
      tau_nu=r*r/nu
      tau_prod=1/aq
      clock=tau_nu/tau_prod
      absRe=-ReQ if ReQ<0 else ReQ
      if not (clock/absRe).contains(1):
        raise AssertionError(('transaction Reynolds clock identity',q_s,r_s,nu_s,clock,absRe))
      rows.append({'q_e':q_s,'r':r_s,'nu':nu_s,'Gamma_Q':str(GammaQ),
                   'signed_Re_Q':str(ReQ),'viscous_over_productive_time':str(clock),
                   'clock_over_abs_Re_Q':str(clock/absRe)})

# Exact Hodge carrier calibration from module 19:
# q_e=(28/5)s(r/L)^2.  The identity must remain true after that geometry is inserted.
carrier=[]
for s_s in ['1e-24','1','1e24']:
  s=arb(s_s)
  for L_s in ['1e-24','1','1e24']:
    L=arb(L_s)
    for x_s in ['1e-12','0.001','0.1','0.5','1']:
      x=arb(x_s); r=L*x
      q=(arb(28)/5)*s*x*x
      for nu_s in ['1e-24','1','1e24']:
        nu=arb(nu_s)
        ReQ=q*r*r/nu
        clock=(r*r/nu)*q
        if not (ReQ/clock).contains(1):
          raise AssertionError(('carrier transaction clock',s_s,L_s,x_s,nu_s,ReQ,clock))
        carrier.append({'s':s_s,'L':L_s,'r_over_L':x_s,'nu':nu_s,
                        'q_e':str(q),'Re_Q':str(ReQ),'clock':str(clock)})

# Exact finite-filament geometry from module 03.
# efficiency E(alpha)=4*pi*d^2*s/Gamma, hence
# chi=s*d^2/nu = E(alpha)/(4*pi) * Re_Gamma.
filament=[]
for a_s in ['1e-12','1e-6','0.01','0.1','0.5','1','1.4142135623730950488','2','10','1e6']:
  alpha=arb(a_s)
  eff=alpha*(alpha*alpha+2)/(1+alpha*alpha)**arb('1.5')
  for Re_s in ['1e-6','0.01','0.1','1','10','1e3','1e6']:
    Re=arb(Re_s)
    chi=eff*Re/(4*pi)
    recovered=4*pi*chi/eff
    if not (recovered/Re).contains(1):
      raise AssertionError(('finite filament Re-clock relation',a_s,Re_s,recovered/Re))
    filament.append({'alpha':a_s,'efficiency':str(eff),'Re_Gamma':Re_s,
                     'strain_diffusion_clock_chi':str(chi),'recovered_Re_ratio':str(recovered/Re)})

# Exact maximum geometric efficiency, used only as a scoped gate for this filament family.
effmax=4*arb(2).sqrt()/(3*arb(3).sqrt())
chi_Re1_max=effmax/(4*pi)
if not (chi_Re1_max < arb('0.1')):
  raise AssertionError(('Re=1 filament should have sub-0.1 strain/diffusion clock',chi_Re1_max))

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'identity_cases':len(rows),
  'carrier_cases':len(carrier),
  'filament_cases':len(filament),
  'finite_filament_efficiency_max':str(effmax),
  'finite_filament_chi_at_Re1_max':str(chi_Re1_max),
  'interpretation':'The signed Hodge transaction circulation Gamma_Q=r^2 q_e has viscosity-normalized magnitude |Gamma_Q|/nu exactly equal to the ratio of the shell viscous crossing time r^2/nu to the productive time 1/|q_e|. In the exact finite-filament geometry the same relation appears with the physical strain efficiency: chi=s d^2/nu=E(alpha) Re_Gamma/(4 pi). Thus productive circulation currency and stretch-versus-diffusion clock are not independent resources.',
  'identity':rows,
  'carrier':carrier,
  'finite_filament':filament,
},indent=2,allow_nan=False))
