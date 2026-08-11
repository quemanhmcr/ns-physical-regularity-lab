import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
# Write epsilon=rho^5 so the validated source horizon R~epsilon^(4/5) is R=K rho^4 exactly.
# maintenance s=mu nu/epsilon^2; fixed material circulation Gamma=ReG nu.
# validated winding transaction scaling: s=Cw Gamma N/R^2.
# Adversarial saturated packing of N passes of the SAME tube: N a^2=R^2.
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
      nu=arb(nus); Gamma=ReG*nu; s=mu*nu/(eps*eps)
      N=s*R*R/(Cw*Gamma)
      a=R/N.sqrt()
      tau_diff=a*a/nu; tau_maint=1/s
      omega=Gamma/(a*a)
      ell=N*R
      volume=a*a*ell
      Ecollar_scale=Gamma*Gamma*ell
      Zscale=omega*omega*volume
      D_event=nu*Zscale*tau_maint
      checks={
       'winding_law':s/(Cw*Gamma*N/(R*R)),
       'packing':N*a*a/(R*R),
       'a_over_eps':(a/eps)/(Cw*ReG/mu).sqrt(),
       'diffusion_buffer':(tau_diff/tau_maint)/(Cw*ReG),
       'vorticity_match':omega/(s/Cw),
       'line_length':ell/(mu*K**3/(Cw*ReG)*rho**2),
       'volume':volume/(K**3*rho**12),
       'collar_scale':Ecollar_scale/(mu*K**3*ReG*nu*nu/Cw*rho**2),
       'event_dissipation':D_event/(mu*K**3*nu*nu/(Cw*Cw)*rho**2),
      }
      for name,val in checks.items():
       if not val.contains(1): raise AssertionError((name,rs,mus,Ks,Cs,Res,nus,val))
      rows.append({
       'rho_epsilon_power_one_fifth':rs,'epsilon':str(eps),'mu_maintenance_coeff':mus,
       'K_source_horizon_coeff':Ks,'Cw_winding_efficiency':Cs,'circulation_Re_Gamma_over_nu':Res,'nu':nus,
       'source_radius_R':str(R),'maintenance_strain_s':str(s),'same_lineage_winding_factor_N':str(N),
       'packed_tube_radius_a':str(a),'a_over_epsilon':str(a/eps),
       'tube_diffusion_over_maintenance_time':str(tau_diff/tau_maint),
       'tube_vorticity_scale_Gamma_over_a2':str(omega),'tube_vorticity_over_maintenance_strain':str(omega/s),
       'same_lineage_length_inside_source':str(ell),'packed_tube_volume_scale':str(volume),
       'fixed_aspect_collar_energy_scale_Gamma2_length':str(Ecollar_scale),
       'enstrophy_scale':str(Zscale),'viscous_dissipation_over_one_maintenance_time':str(D_event),
       'all_closed_identity_ratios':{k:str(v) for k,v in checks.items()}})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'Use the validated maintenance scaling s=mu nu/epsilon^2 and finite-energy source horizon R=K epsilon^(4/5), writing epsilon=rho^5. Let one material lineage keep fixed circulation Gamma=Re_Gamma nu and supply the source by the validated winding law s=Cw Gamma N/R^2. Then N=(mu K^2/(Cw Re_Gamma)) rho^(-2) diverges as epsilon->0 without changing Gamma. '
  'Give the N passes the largest tube radius allowed by saturated cross-sectional packing, N a^2=R^2. Exact cancellation gives a/epsilon=sqrt(Cw Re_Gamma/mu) and (a^2/nu)/(1/s)=Cw Re_Gamma, independent of epsilon. Thus a fixed high-Re lineage retains a fixed viscous-persistence buffer while its winding multiplicity diverges. Its vorticity scale Gamma/a^2 equals s/Cw, so the same lineage is itself amplified to the current maintenance scale rather than being a weak reusable wire. '
  'The required lineage length inside the source is proportional to rho^2=epsilon^(2/5), the packed tube volume to rho^12, fixed-aspect circulation-collar occupancy Gamma^2 ell to epsilon^(2/5), and viscous enstrophy dissipation accumulated over one maintenance time also to epsilon^(2/5). These scalar costs vanish at the shrinking end for fixed coefficients. '
  'This is an adversarial scaling construction, not an exact Navier-Stokes blow-up solution: saturated clean packing, the winding law coefficient, and continual stretching of the same material tube are assumed. It KILLs only the naive theorem that infinitely many maintenance events or diverging fixed-time winding necessarily consume infinitely many distinct circulation packets, or that packing plus tube diffusion alone forbids reuse. The remaining causal bottleneck is who continually generates the winding/stretch of this already-high-Re lineage while respecting closure and the directed interaction network.'),
 'rows':rows},indent=2,allow_nan=False))
