import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
for n in (3,5,9,17,33,65,129,257):
 N=arb(n)
 cG=2*(2*N+1)/(N*(N+1)*(N-1)*(N-2))
 cClock=2*N*(2*N+1)/((N-1)*(N-2))
 for ls in ('1','1.01','1.05','1.1','1.25','1.5','2'):
  lam=arb(ls)
  qgain=lam**(n-3)
  threshold=cClock*(lam**(n-2))
  Hcrit=N*N*(N+1)/(lam*lam)
  Rcrit=cG*Hcrit*(lam**n)
  if not (Rcrit/threshold).contains(1): raise AssertionError(('threshold',n,ls))
  for Hs in ('1e-30','1','1e30'):
   H=arb(Hs)
   Rg=cG*H*(lam**n)
   clock=H*(lam*lam)/(N*N*(N+1))
   conv=cClock*(lam**(n-2))
   if not (Rg/(conv*clock)).contains(1): raise AssertionError(('identity',n,ls,Hs))
   rows.append({'n_degree':n,'lambda_source_radius_over_core':ls,'H_core_Qeps_eps3_over_nu':Hs,
    'curvature_source_gain_Q_L_over_Q_eps':str(qgain),'source_lobe_circulation_Re':str(Rg),
    'source_angular_clock_over_core_one_wavelength_fold_clock':str(clock),
    'exact_Re_per_clock_conversion':str(conv),'near_contact_conversion_lambda1':str(cClock),
    'conversion_over_near_contact':str(conv/cClock),'H_core_required_for_clock_equal_1':str(Hcrit),
    'Gamma_source_Re_at_clock_equal_1':str(Rcrit),'closed_persistence_threshold_Re':str(threshold),
    'threshold_over_near_contact':str(threshold/cClock)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'Continue the matched harmonic folding source from core epsilon to L=lambda epsilon. Homogeneity gives Q_L/Q_epsilon=lambda^(n-3), so Gamma_source/nu=c_G H_epsilon lambda^n with H_epsilon=Q_epsilon epsilon^3/nu. '
  'The source tangential angular clock grows only as lambda^2. Eliminating H_epsilon gives Gamma_source/nu=[2n(2n+1)/((n-1)(n-2))] lambda^(n-2) times the source-angular/core-fold clock ratio. '
  'Thus survival for one core one-wavelength folding time forces Gamma_source/nu>=C_n lambda^(n-2), C_n->4. For lambda>=1 the minimum threshold is near contact; remote placement worsens it exponentially in degree. '
  'This is a single-degree matched-source persistence gate, not a universal 3D theorem. Continuous nonlinear regeneration, mixtures, and repeated low-degree folding remain explicit escapes.'),
 'rows':rows},indent=2,allow_nan=False))
