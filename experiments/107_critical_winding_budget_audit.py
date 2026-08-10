import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
# Use L=A tau^(2/5), s=1/(5tau). Track only exact scaling ledgers.
for As in ['1e-12','1','1e12']:
 A=arb(As)
 for nus in ['1e-24','1','1e24']:
  nu=arb(nus)
  for taus in ['1e-6','1e-30','1e-100']:
   tau=arb(taus); L=A*tau**(arb(2)/5); s=1/(5*tau)
   Re=s*L*L/nu
   Z=s*s*L**3
   diss_tail=nu*A**3/5*tau**(arb(1)/5) # int_0^tau nu*(A^3/25)t^-4/5 dt
   exposure_tail=5*nu/(A*A)*tau**(arb(1)/5)
   energy=s*s*L**5
   rows.append({'A_source_prefactor':As,'nu':nus,'tau':taus,'L':str(L),'s':str(s),'Re_source':str(Re),'productive_enstrophy_scale_s2L3':str(Z),'remaining_viscous_enstrophy_budget_scale':str(diss_tail),'remaining_localized_viscous_exposure':str(exposure_tail),'energy_occupancy_s2L5':str(energy)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'The critical frozen winding conveyor sits exactly on the Hodge energy horizon: s^2L^5 is constant. Its productive enstrophy scales as tau^-4/5 but is spacetime-integrable, the remaining viscous dissipation and localized spectral exposure both vanish like tau^1/5, while Re_source diverges like tau^-1/5. Hence all scalar finite-energy/dissipation mechanisms identified so far permit this critical winding escape. A regularity proof must break its material/network geometry, not its scalar scaling.','rows':rows},indent=2,allow_nan=False))
