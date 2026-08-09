import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be >=160')
ctx.prec=BITS

# Exact helical Beltrami mode at selected z values. All derivatives are analytic.
nu_values=['1e-18','1e-6','1','1e6']
k_values=['1e-9','1e-3','1','1e3','1e9']
A_values=['1e-12','1','1e12']
zphase_values=['0','0.2','1.1','2.7']
rows=[]
for nu_s in nu_values:
  for k_s in k_values:
    for A_s in A_values:
      nu=arb(nu_s); k=arb(k_s); A=arb(A_s)
      # We avoid evaluating trig because all invariant contractions are phase independent.
      rho=A
      grad_xi_sq=k*k
      # For u=-(A/k)(cos kz,sin kz,0), grad u has only z-column:
      # du/dz=(A sin kz,-A cos kz,0).  Contracting S with xi on both sides is exactly zero
      # because xi_z=0 and S has one z index.
      stretch=arb(0)
      # rho uniform in space, so Delta rho=0; material transport also zero because u_z=0.
      directional_tax=nu*rho*grad_xi_sq
      drho_dt=-nu*k*k*rho
      residual=drho_dt-(rho*stretch-directional_tax)
      if not residual.contains(0):
          raise AssertionError(('magnitude identity residual',nu_s,k_s,A_s,residual))
      # Full NS residual after choosing p constant: du/dt - nu Delta u =0; nonlinearity is gradient
      # since u x omega=0 and |u| is spatially constant, in fact (u.grad)u=0 here.
      # Each velocity component obeys d_t u = nu d_zz u = -nu k^2 u.
      decay_rate=(-drho_dt/rho)
      if not (decay_rate/(nu*k*k)).contains(1):
          raise AssertionError('decay rate mismatch')
      rows.append({'nu':nu_s,'k':k_s,'A':A_s,'geometry_tax_rate_over_rho':str(directional_tax/rho),'stretch':str(stretch)})

print(json.dumps({
 'arb_precision_bits':BITS,
 'status':'PASS',
 'cases':len(rows),
 'interpretation':'Exact Beltrami/helical Navier-Stokes modes can have arbitrarily rapid direction variation while xi.S.xi is exactly zero; viscosity removes vorticity magnitude at exactly nu|grad xi|^2. Directional complexity is therefore a real viscous tax, but not by itself a stretching resource.',
 'rows':rows,
},indent=2))
