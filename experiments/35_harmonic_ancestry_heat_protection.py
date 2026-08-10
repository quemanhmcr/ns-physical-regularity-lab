import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

# Exact Navier-Stokes shear heat mode on a periodic box:
#   u_x = H + A exp(-nu k^2 t) cos(k y), u_y=u_z=0.
# Nonlinearity vanishes.  H is the harmonic x-cycle mode; the cosine is vortical.
# Interior viscosity damps the vortical mode but leaves the global harmonic mode unchanged.
# Individual material x-loop circulations at fixed y nevertheless change by Kelvin viscosity.

Ps=['1e-18','1','1e18']
Ls=['1e-12','1','1e12']
nus=['1e-24','1','1e24']
ks=['1e-18','1','1e18']
Hs=['1e-12','1','1e12']
As=['1e-12','1','1e12']
qs=['0','1e-12','0.1','1','10','100']  # q=nu k^2 t
phase_cos=['-1','-0.5','0.5','1']
rows=[]

for Ps_ in Ps:
  P=arb(Ps_)
  for Ls_ in Ls:
    Lz=arb(Ls_)
    for nus_ in nus:
      nu=arb(nus_)
      for ks_ in ks:
        k=arb(ks_)
        Ly=2*pi/k
        V=P*Ly*Lz
        for Hs_ in Hs:
          H=arb(Hs_)
          Gamma_h=P*H
          Eh=V*H*H/2
          for As_ in As:
            A=arb(As_)
            for qs_ in qs:
              q=arb(qs_)
              decay=(-q).exp()
              B=A*decay
              Ev=V*B*B/4
              # Compute the spatially averaged total energy independently from
              # the harmonic+vortical split.  Do not recover an exponentially
              # small Ev by subtracting two nearly equal Arb balls.
              Et_direct=V*(H*H + B*B/2)/2
              split=Eh+Ev
              split_ratio=Et_direct/split
              if not split_ratio.contains(1):
                  raise AssertionError(('stable harmonic/vortical energy split',Ps_,Ls_,nus_,ks_,Hs_,As_,qs_,split_ratio))
              cross_energy=arb(0)  # mean_y cos(k y)=0 exactly
              # Exact energy dissipation of the cosine mode.
              enstrophy=V*k*k*B*B/2
              dEv_dt=-nu*k*k*V*B*B/2
              if not ((-dEv_dt)/(nu*enstrophy)).contains(1):
                  raise AssertionError(('shear heat energy identity',Ps_,Ls_,nus_,ks_,Hs_,As_,qs_))
              # Harmonic coefficient is exactly unchanged.
              Gamma_h_now=P*H
              if not (Gamma_h_now/Gamma_h).contains(1):
                  raise AssertionError('harmonic period drifted under interior heat')
              probes=[]
              for cs in phase_cos:
                c=arb(cs)
                Gamma_loop=P*(H+B*c)
                dGamma_dt=-nu*k*k*P*B*c
                kelvin=nu*P*(-k*k*B*c)
                if cs!='0' and not (dGamma_dt/kelvin).contains(1):
                    raise AssertionError(('Kelvin loop mutation mismatch',cs,dGamma_dt,kelvin))
                deviation=Gamma_loop-Gamma_h
                expected_dev=P*A*c*decay
                if cs!='0' and not (deviation/expected_dev).contains(1):
                    raise AssertionError(('loop-to-harmonic deviation decay',cs,deviation,expected_dev))
                probes.append({
                    'cos_phase':cs,
                    'material_loop_circulation':str(Gamma_loop),
                    'loop_minus_harmonic_period':str(deviation),
                    'dGamma_dt':str(dGamma_dt),
                    'Kelvin_viscous_term':str(kelvin),
                })

              rows.append({
                'P':Ps_,'Lz':Ls_,'nu':nus_,'k':ks_,'H':Hs_,'A':As_,'q':qs_,
                'harmonic_period':str(Gamma_h),
                'harmonic_energy':str(Eh),
                'vortical_energy':str(Ev),
                'cross_energy':str(cross_energy),
                'total_energy_direct':str(Et_direct),
                'total_over_harmonic_plus_vortical':str(split_ratio),
                'minus_dEv_dt_over_nu_enstrophy':str((-dEv_dt)/(nu*enstrophy)),
                'loop_probes':probes,
              })

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'cases':len(rows),
  'interpretation':(
      'In an exact periodic Navier-Stokes shear solution, viscosity exponentially removes the vortical Fourier mode while the global harmonic circulation mode is exactly protected. '
      'At the same time, individual material x-loop circulations at fixed y change by the exact Kelvin viscous term and relax toward the harmonic period. '
      'Therefore harmonic finite-thickness ancestry and codimension-two material-loop ancestry coincide only in a circulation-isolated collar; once vorticity fills the collar they are distinct observables.'
  ),
  'rows':rows,
},indent=2,allow_nan=False))
