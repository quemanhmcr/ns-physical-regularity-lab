import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
one=arb(1); two=arb(2); pi=arb.pi(); rt2=arb(2).sqrt()
# Eigenvalues of the validated stationary amplifier strain S_*/sigma.
c_plus=(3*rt2-one)/2          # largest stretching eigenvalue
c_minus=(one+3*rt2)/2         # magnitude of most compressive eigenvalue
# Remaining eigenvalue is 1; trace is zero.
if not (c_plus+one-c_minus).contains(0): raise AssertionError('trace eigenvalue identity')
if not (c_plus*c_minus-arb('4.25')).contains(0): raise AssertionError('eigen product identity')

# Along the most compressive eigenvector, k grows as exp(c_minus sigma t).
# Choose Fourier vorticity polarization along the most stretching eigenvector (orthogonal because S_* is symmetric).
# Then exact characteristic amplitude is exp[c_plus sigma t - nu M (exp(2 c_minus sigma t)-1)/(2 c_minus sigma)].
def log_gain_one_initial_viscous_clock(x):
    # x=sigma/(nu M), t=1/(nu M)
    return c_plus*x - ((2*c_minus*x).exp()-one)/(2*c_minus*x)

# Closed strongest-cascade characteristic sampled only as an autopsy.  No scalar x is promoted as a universal optimum.
x_samples=[]
for xs in ['0.1','0.25','0.5','0.75','1','1.5','2','3']:
    x=arb(xs)
    x_samples.append({'sigma_over_nuM':xs,'one_initial_viscous_clock_log_gain':str(log_gain_one_initial_viscous_clock(x))})
x_inst=one/c_plus
inst_one_clock_gain=log_gain_one_initial_viscous_clock(x_inst)

rows=[]
E0=one; nu=one
for Ms in ['1e2','1e4','1e8','1e16','1e32','1e64','1e128']:
    M=arb(Ms)
    sigma_inst=nu*M/c_plus
    # Validated stationary-amplifier harmonic occupancy: E_h=(7 pi/5) sigma^2 R^5.
    R_inst=(arb(5)*E0/(arb(7)*pi*sigma_inst*sigma_inst))**(one/5)
    eps=one/M.sqrt()
    rows.append({
        'M_k_squared':Ms,
        'collar_length_epsilon_M_minus_half':str(eps),
        'instantaneous_nondecay_min_sigma_nu1':str(sigma_inst),
        'instantaneous_horizon_identity_nu_M_over_sigma':str(nu*M/sigma_inst),
        'instantaneous_horizon_sigma_over_nuM':str(x_inst),
        'strongest_cascade_log_gain_after_one_initial_viscous_clock_when_entering_at_instantaneous_horizon':str(inst_one_clock_gain),
        'finite_energy_stationary_strain_radius_ceiling_instantaneous':str(R_inst),
        'R_inst_times_M_power_2_over_5':str(R_inst*(M**(arb(2)/5))),
        'R_inst_over_collar_epsilon':str(R_inst/eps),
    })

# Scaling gate independent of the chosen stationary eigenframe:
# for omega_epsilon=A_epsilon Omega(x/epsilon), E_S omega scales like A_epsilon,
# while nu Delta omega scales like nu A_epsilon epsilon^-2.
scaling=[]
for es in ['1','1e-2','1e-4','1e-8','1e-16','1e-32']:
    e=arb(es)
    scaling.append({'epsilon':es,'affine_to_viscous_operator_scale_ratio_per_unit_sigma_over_nu':str(e*e),'required_sigma_over_nu_for_order_one_balance':str(one/(e*e))})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS',
 'stationary_strain_eigenvalues_over_sigma':[str(c_plus),'1',str(-c_minus)],
 'instantaneous_vorticity_stretch_coefficient_c_plus':str(c_plus),
 'compressive_wavevector_growth_coefficient_c_minus':str(c_minus),
 'instantaneous_horizon_sigma_over_nuM':str(x_inst),
 'strongest_cascade_characteristic_autopsy':x_samples,
 'interpretation':(
   'For the validated stationary productive amplifier S_*, the symmetric strain eigenvalues are c_plus sigma, sigma and -c_minus sigma with c_plus=(3sqrt2-1)/2 and c_minus=(1+3sqrt2)/2. '
   'In an affine Navier-Stokes calibration, a Fourier vorticity packet polarized along the strongest stretching eigenvector and with wavevector along the strongest compressive eigenvector obeys exactly k(t)=k0 exp(c_minus sigma t) and log |omega_hat(t)/omega_hat(0)|=c_plus sigma t-nu k0^2[exp(2c_minus sigma t)-1]/(2c_minus sigma). '
   'The instantaneous nondecay gate is nu k^2 <= c_plus sigma.  Thus a remote collar at k^2=M cannot even avoid instantaneous viscous loss under a fixed reusable affine strain as M tends to infinity; the required strain rate itself must scale like nu M.  The exact strongest-cascade characteristic is reported only as an autopsy because increasing sigma also accelerates wavevector growth and is not a one-parameter universal optimization. '
   'The generic scaling identity E_S[Omega(x/epsilon)]~sigma Omega while nu Delta[Omega(x/epsilon)]~nu epsilon^-2 Omega gives the same gate sigma epsilon^2/nu=O(1). '
   'Combining sigma~nu M with the already validated stationary-amplifier Hodge occupancy E_h=(7pi/5)sigma^2 R^5 yields R_source=O(M^-2/5)=O(epsilon^(4/5)). '
   'The Hodge exponent 4/5 therefore reappears from the actual requirement to replenish a remote viscous collar.  This is not a contradiction: R_source/epsilon grows like epsilon^-1/5.  It says that the cheap static collar escape forces an inward cascade of increasingly strong strain sources rather than being maintainable by one fixed large-scale affine reservoir.'),
 'rows':rows,'generic_affine_scaling_rows':scaling
},indent=2,allow_nan=False))
