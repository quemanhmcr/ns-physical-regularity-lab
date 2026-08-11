import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
one=arb(1); rows=[]
# Matched degree-n harmonic source sphere.  Odd n=m+1 corresponds to modules 184-185.
# phi_in=-(a/n)r^nY_n, phi_out=a L^(2n+1)/(n+1) r^(-n-1)Y_n.
# Normal velocity matches at r=L.  Equatorial jump is
# [u_theta]=-a(2n+1)/(n+1)L^(n-1)sin(n theta).
for n in (3,5,9,17,33,65,129,257):
 N=arb(n)
 cG=2*(2*N+1)/(N*(N+1)*(N-1)*(N-2))
 cClock=2*N*(2*N+1)/((N-1)*(N-2))
 radial_in=-one
 radial_out=-(N+1)*(one/(N+1))
 if not (radial_in-radial_out).contains(0): raise AssertionError(('normal match',n,radial_in,radial_out))
 # Direct Stokes reconstruction for a=L=1.
 jump=(2*N+1)/(N+1)
 int_lobe=2/N
 Gdirect=jump*int_lobe
 QL=(N-1)*(N-2)
 Gclosed=cG*QL
 if not (Gdirect/Gclosed).contains(1): raise AssertionError(('Gamma',n,Gdirect,Gclosed))
 for Hs in ('1e-60','1e-12','1','1e12','1e60'):
  H=arb(Hs) # H=Q_L L^3/nu
  Rg=cG*H
  clock=H/(N*N*(N+1)) # tau_ang/tau_fold
  if not (Rg/(cClock*clock)).contains(1): raise AssertionError(('clock',n,Hs))
  rows.append({
   'n_spherical_harmonic_degree':n,'m_velocity_axis_power':n-1,
   'H_curvature_Re_Q_L3_over_nu':Hs,
   'normal_match_phi_out_coefficient_over_a_L2n1':str(one/(N+1)),
   'a1_L1_radial_velocity_inner':str(radial_in),'a1_L1_radial_velocity_outer':str(radial_out),
   'normal_velocity_match_residual':str(radial_in-radial_out),
   'one_positive_lobe_integral_abs_sin_n_theta':str(int_lobe),
   'Gamma_lobe_over_Q_L3':str(cG),'circulation_Re_Gamma_lobe_over_nu':str(Rg),
   'tau_angular_diffusion_over_tau_one_wavelength_fold':str(clock),
   'exact_GammaRe_over_clock_factor':str(cClock),
   'factor_over_asymptotic_4':str(cClock/4),
   'direct_Gamma_a1_L1':str(Gdirect),'closed_Gamma_a1_L1':str(Gclosed),
   'direct_over_closed_Gamma':str(Gdirect/Gclosed)})
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'Match the degree-n interior harmonic folding potential phi_in=-(a/n)r^nY_n to the unique decaying exterior harmonic potential with identical normal velocity at r=L. The piecewise velocity is distributionally incompressible and its source vorticity is a tangential degree-n vortex sheet. '
  'One positive equatorial jump lobe has actual Stokes circulation Gamma_lobe=2|a|(2n+1)L^n/[n(n+1)]. With boundary curvature source Q_L=|a|(n-1)(n-2)L^(n-3), Gamma_lobe/(Q_L L^3)=2(2n+1)/[n(n+1)(n-1)(n-2)]~4/n^3. Thus static circulation per lobe can vanish at high degree. '
  'The degree-n tangential source has intrinsic angular viscous clock tau_ang=L^2/[nu n(n+1)]. Define tau_fold=(n/L)/Q_L, the time to build curvature equal to the inverse angular wavelength. Exactly Gamma_lobe/nu=[2n(2n+1)/((n-1)(n-2))](tau_ang/tau_fold), with factor tending to 4. '
  'Therefore high degree does not require a fixed order-nu static circulation packet, but Gamma_lobe/nu->0 makes angular diffusion faster than one-wavelength folding in this thin-shell calibration. The vortex sheet is an ideal limit; smooth shell regeneration, degree mixtures, and repeated low-degree folding remain escapes.'),
 'rows':rows},indent=2,allow_nan=False))
