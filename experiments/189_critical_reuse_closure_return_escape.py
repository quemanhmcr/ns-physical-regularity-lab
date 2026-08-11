import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); one=arb(1); rows=[]; skipped_nonmultipass=0
# epsilon=rho^5; source R=K rho^4; maintenance s=mu nu/rho^10.
# Fixed lineage Gamma=ReG nu and winding law s=Cw Gamma N/R^2.
# Saturated passage packing gives a=R/sqrt(N)=sqrt(Cw ReG/mu) rho^5.
# Connect productive passes by either:
#   near hairpins Rc=chi a,
#   remote returns Rc=eta R.
# N semicircular returns have total bend length Lb=N*pi*Rc and
# curvature-square inventory K2=N*pi/Rc.
# A fixed-aspect circulation isolation collar b=beta a gives the exact tubular
# Cauchy floor Ecoll >= theta^2 Gamma^2 N Rc/4 * (1-b/Rc) log(beta),
# whenever b<Rc.  Directional viscous curvature rate is nu Gamma K2.
# We audit its one-maintenance-time integral and the length-averaged curvature
# diffusion exposure nu tau K2/Lb=nu tau/Rc^2.
for rs in ('0.1','1e-3','1e-6'):
 rho=arb(rs); eps=rho**5
 for mus in ('0.3','1','7'):
  mu=arb(mus)
  for Ks in ('0.5','2'):
   K=arb(Ks); R=K*rho**4
   for Cs in ('0.2','1.3'):
    Cw=arb(Cs)
    for Res in ('10','1e3','1e6'):
     ReG=arb(Res)
     for nus in ('1e-9','1','1e9'):
      nu=arb(nus); Gamma=ReG*nu
      s=mu*nu/(eps*eps); tau=1/s
      N=s*R*R/(Cw*Gamma)
      # A return geometry is meaningful only in the multi-pass branch.
      # The original run 31477455263 incorrectly included coarse cases N<1.
      if not (N-one).lower()>0:
       skipped_nonmultipass += 1
       continue
      a=R/N.sqrt()
      # Two geometric closure extremes. beta<chi for the near tubular Jacobian.
      for kind,fac_s,beta_s,theta_s in (
       ('near_hairpin','3','1.5','0.7'),
       ('remote_return','2','1.5','0.7')):
       fac=arb(fac_s); beta=arb(beta_s); theta=arb(theta_s)
       Rc=fac*a if kind=='near_hairpin' else fac*R
       b=beta*a
       if not (Rc-b).lower()>0: raise AssertionError(('invalid tubular collar',kind,rs,mus,Ks,Cs,Res,nus,Rc,b))
       Lbend=N*pi*Rc
       K2=N*pi/Rc
       Efloor=theta*theta*Gamma*Gamma*N*Rc/4*(1-b/Rc)*beta.log()
       tax_rate=nu*Gamma*K2
       tax_event=tax_rate*tau
       avg_exposure=nu*tau*K2/Lbend
       # Structural identities independent of the particular asymptotic branch.
       checks={
        'semicircle_length_curvature_product': Lbend*K2/(N*N*pi*pi),
        'average_curvature_identity': avg_exposure/(nu*tau/(Rc*Rc)),
        'collar_floor_reconstruction': Efloor/(theta*theta*Gamma*Gamma*N*Rc/4*(1-b/Rc)*beta.log()),
        'tax_reconstruction': tax_event/(nu*Gamma*N*pi/Rc*tau),
       }
       if kind=='near_hairpin':
        # Exact critical cancellation: exposure is scale independent and high-Re small.
        closed_exp=1/(fac*fac*Cw*ReG)
        closed_L=pi*fac*K*K*(mu/(Cw*ReG)).sqrt()*rho**3
        closed_tax=pi*nu*K*K/(fac*Cw)*(mu/(Cw*ReG)).sqrt()*rho**3
        checks['near_exposure']=avg_exposure/closed_exp
        checks['near_length_rho3']=Lbend/closed_L
        checks['near_event_tax_rho3']=tax_event/closed_tax
       else:
        # Remote returns have even smaller average exposure as rho->0.
        closed_exp=rho**2/(mu*fac*fac*K*K)
        closed_L=pi*fac*mu*K**3/(Cw*ReG)*rho**2
        closed_tax=pi*nu*K/(fac*Cw)*rho**4
        checks['remote_exposure_rho2']=avg_exposure/closed_exp
        checks['remote_length_rho2']=Lbend/closed_L
        checks['remote_event_tax_rho4']=tax_event/closed_tax
       for name,val in checks.items():
        if not val.contains(1): raise AssertionError((name,kind,rs,mus,Ks,Cs,Res,nus,val))
       rows.append({
        'rho_epsilon_power_one_fifth':rs,'epsilon':str(eps),'closure_kind':kind,
        'mu_maintenance_coeff':mus,'K_source_horizon_coeff':Ks,'Cw_winding_efficiency':Cs,
        'circulation_Re_Gamma_over_nu':Res,'nu':nus,'return_radius_factor':fac_s,
        'collar_outer_over_tube_radius_beta':beta_s,'collar_circulation_fraction_theta':theta_s,
        'source_radius_R':str(R),'winding_factor_N':str(N),'packed_tube_radius_a':str(a),
        'return_curvature_radius_Rc':str(Rc),'return_total_bend_length':str(Lbend),
        'return_curvature_square_inventory':str(K2),'fixed_aspect_collar_energy_floor':str(Efloor),
        'directional_curvature_viscous_rate':str(tax_rate),
        'directional_curvature_tax_over_one_maintenance_time':str(tax_event),
        'length_averaged_curvature_diffusion_exposure':str(avg_exposure),
        'all_closed_identity_ratios':{k:str(v) for k,v in checks.items()}})
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'skipped_nonmultipass_parameter_points':skipped_nonmultipass,
 'interpretation':(
  'Close the module-188 same-lineage productive passes by N semicircular returns and audit the exact curvature-square inventory together with the tubular circulation-collar floor already derived in the closure microscope. Two adversarial limits are tested: near hairpins of radius chi a and remote returns of radius eta R. '
  'For near hairpins, total bend length scales as epsilon^(3/5), the curvature-viscous tax integrated over one maintenance time also scales as epsilon^(3/5), while the length-averaged curvature diffusion exposure is exactly 1/(chi^2 Cw Re_Gamma), independent of epsilon. Thus a fixed high-Re lineage can make the tight-return branch persist on each shrinking maintenance time even though total curvature-square diverges. '
  'For remote returns, total bend length scales as epsilon^(2/5), integrated curvature tax as epsilon^(4/5), and the average curvature exposure as epsilon^(2/5). Remote closure is even less diffusive on the shrinking maintenance clock. The fixed-aspect circulation-collar floor is positive and audited with its tubular Jacobian, but its scale follows the shrinking return length rather than producing a non-summable packet toll. '
  'This KILLs the hope that closure curvature plus a local circulation collar automatically repairs the packet-counting theorem at the critical source exponents. It is still a geometric scaling construction: it does not construct a globally embedded non-self-intersecting tube, control collective Biot-Savart induction, or show that the return geometry is dynamically maintained by Navier-Stokes. Those actor-of-actor dynamics remain the bottleneck.'),
 'rows':rows},indent=2,allow_nan=False))
