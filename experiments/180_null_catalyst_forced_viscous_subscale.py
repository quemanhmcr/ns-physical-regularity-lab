import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
# Critical branch: s=tau^-1, epsilon=tau^alpha, beta=(1+4alpha)/7.
# Null catalyst gradient amplitude g ~ sqrt(s nu)/epsilon^2 ~ sqrt(nu) tau^(-7beta/2).
# If B=g B0 keeps its spectral ray, the jet law requires g_dot=nu c3, with c3 a grad Delta omega coefficient.
# Define curvature length ell^2=g/|c3|. Exact required ell^2=2 nu tau/(7 beta).
for alphas in ('0.4','0.45','0.49'):
 alpha=arb(alphas);beta=(1+4*alpha)/7
 for nus in ('1e-12','1','1e12'):
  nu=arb(nus)
  for taus in ('1e-6','1e-30','1e-100'):
   tau=arb(taus);eps=tau**alpha;strain=1/tau;Rs=strain*eps*eps/nu
   g=(strain*nu).sqrt()/(eps*eps)
   growth=(7*beta/2)*g/tau
   c3=growth/nu
   ell2=g/c3;ell=ell2.sqrt()
   pred2=2*nu*tau/(7*beta)
   if not (ell2/pred2).contains(1):raise AssertionError(('ell2',alphas,nus,taus,ell2,pred2))
   ratio=ell/eps;predratio=(arb(2)/(7*beta*Rs)).sqrt()
   if not (ratio/predratio).contains(1):raise AssertionError(('ratio',alphas,nus,taus,ratio,predratio))
   rows.append({'alpha':alphas,'beta':str(beta),'nu':nus,'tau':taus,'epsilon':str(eps),'source_Re_s_epsilon2_over_nu':str(Rs),'required_null_gradient_amplitude_scale_g':str(g),'required_material_growth_rate_gdot':str(growth),'required_aligned_grad_Delta_omega_coefficient':str(c3),'forced_viscous_curvature_length_ell':str(ell),'ell_squared':str(ell2),'closed_ell2_2nutau_over_7beta':str(pred2),'ell_over_core_epsilon':str(ratio),'closed_ell_over_epsilon_sqrt_2_over_7betaRs':str(predratio)})
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'Apply the vorticity-zero gradient spectral law to the critical linear-null catalyst.  Its required gradient strength grows as g~sqrt(nu) tau^(-7 beta/2).  Euler commutators cannot change the spectral amplitude, so on a fixed material vorticity-zero center the growth must be supplied by nu grad Delta omega. '
  'If ell is the physical curvature length defined by |grad Delta omega|~g/ell^2 along the required spectral ray, exact balance gives ell^2=2 nu tau/(7 beta).  Relative to the productive core epsilon, ell/epsilon=sqrt[2/(7 beta R_s)], where R_s=s epsilon^2/nu is the source Reynolds. '
  'Therefore in the high-Re branch R_s->infinity, any viscosity-mediated amplification of the null-catalyst eigenvalues requires a new length ell<<epsilon, asymptotically the diffusive time-to-singularity scale sqrt(nu tau) up to an order-one factor. '
  'This routes the null-catalyst escape back into an unavoidable inward viscous subscale.  It is conditional on following the same material vorticity-zero catalyst; switching/recruiting fresh high-gradient zeros remains an explicit escape and must be attacked separately.'),
 'rows':rows
},indent=2,allow_nan=False))
