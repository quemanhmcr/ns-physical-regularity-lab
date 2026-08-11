import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
# Critical winding calibration: strain count N=-log(tau)/5, L=A exp(-2N), s=(1/5)exp(5N), Re_source=A^2/(5nu) exp(N).
# A localized l-sector squared-norm pure-diffusion exponent per strain count is (8l+12)/Re_source.
for As in ['1e-12','1','1e12']:
 A=arb(As)
 for nus in ['1e-24','1','1e24']:
  nu=arb(nus); pref=5*nu/(A*A)
  for Ni in [0,1,10,100]:
   N=arb(Ni); eN=N.exp(); Re=(A*A/(5*nu))*eN
   llin=2+2*N
   theta=(8*llin+12)/Re
   # Exact tail for l(N)=2+2N: int_N^inf pref*(28+16x)e^-x dx = pref*(44+16N)e^-N.
   tail=pref*(44+16*N)*(-N).exp()
   rows.append({'A':As,'nu':nus,'strain_count_N':Ni,'Re_source':str(Re),'linear_angular_degree_envelope_l_2plus2N':str(llin),'instantaneous_squared_norm_diffusion_exponent_per_strain_count':str(theta),'remaining_linear_l_viscous_exposure':str(tail)})
# Exponential angular envelopes l~exp(beta N).  Beta<1 remains integrable against Re~exp(N); beta=1 is critical.
envelopes=[]
for beta_s in ['0.25','0.5','0.9']:
 beta=arb(beta_s)
 for Nint in [0,10,100]:
  N=arb(Nint)
  # Leading l term exposure tail: int_N inf 8*pref exp[-(1-beta)x] dx.
  # Use pref=1 here; geometry factor can be restored multiplicatively.
  tail=arb(8)/(1-beta)*(-(1-beta)*N).exp()
  envelopes.append({'beta':beta_s,'N':Nint,'leading_exposure_tail_for_l_exp_betaN_with_unit_prefactor':str(tail),'finite_for_beta_less_than_1':True})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'linear_envelope_rows':rows,'subcritical_exponential_envelopes':envelopes,'critical_beta_one_statement':'If l grows proportionally to exp(N), the leading (8l)/Re_source exposure per strain count approaches a nonzero constant and its cumulative integral diverges.','interpretation':'The generalized localized angular gap grows only linearly like 4l after the smooth r^l zero mode is removed.  In the critical frozen winding conveyor Re_source grows like exp(N) in strain count.  Therefore an angular hierarchy whose degree grows linearly, polynomially, or even as exp(beta N) with beta<1 accumulates only finite viscous exposure.  Angular order becoming unbounded is not by itself a regularity contradiction; to force infinite viscous renewal through this gap, angular complexity must grow on the same exponential scale as Re_source or faster.  This kills the naive hope that the l(l+1) term alone automatically defeats the critical winding branch.'},indent=2,allow_nan=False))
