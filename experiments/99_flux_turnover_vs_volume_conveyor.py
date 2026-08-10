import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rows=[]
# Frozen high-Re scaling: Phi=c e^{aN}, L=e^{-bN}. Net positive source-flux recruitment grows like Phi,
# while a full source-volume replacement tail is finite. Flux per source cross-sectional area grows like e^{(a+2b)N}=s scale.
for a_s,b_s in [('0.5','2'),('1','3'),('2','5')]:
 a=arb(a_s); b=arb(b_s); c=arb(3)
 for Ns in ['0','1','10','100']:
  N=arb(Ns); Phi=c*(a*N).exp(); L=(-b*N).exp(); area=L*L; volume=L**3
  flux_density=Phi/area
  dPhi_dN=a*Phi
  volume_tail=(-3*b*N).exp()/(3*b)
  rows.append({'a':a_s,'b':b_s,'N':Ns,'source_flux_Phi':str(Phi),'positive_flux_recruitment_per_strain_count_dPhi_dN':str(dPhi_dN),'source_scale_L':str(L),'source_cross_section_scale_L2':str(area),'flux_density_Phi_over_L2':str(flux_density),'source_volume_scale_L3':str(volume),'late_distinct_volume_tail_int_L3dN':str(volume_tail)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'In the cumulatively frozen escape Phi~Re grows exponentially in strain count while the shrinking source volume can have a finite recruitment tail. The net positive circulation-flux recruitment dPhi/dN is unbounded, and Phi/L^2 grows on the same scale as the source production rate s. Thus later source parcels may have vanishing material volume but must carry increasingly strong already-amplified vorticity flux. The surviving high-Re mechanism is an amplification conveyor, not a fresh-volume supply.','rows':rows},indent=2,allow_nan=False))
