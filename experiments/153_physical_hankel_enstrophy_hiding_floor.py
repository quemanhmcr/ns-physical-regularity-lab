import json, os, math
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

rows=[]
for l in (1,2,4,6,8,12):
    d=2*l+3
    nu=arb(d)/2-1
    # Self-inverse/unitary effective-dimensional radial Hankel transform:
    # a(r)=c_d int A(k) Phi_d(kr) k^(d-1) dk,
    # c_d=[2^(d/2-1) Gamma(d/2)]^-1.
    cd=one=arb(1)
    cd=one/((arb(2)**nu)*(arb(d)/2).gamma())
    gaussian_radial_integral=(arb(2)**(arb(d)/2-1))*(arb(d)/2).gamma()
    gaussian_normalization=cd*gaussian_radial_integral
    if not gaussian_normalization.contains(1):
        raise AssertionError(('Hankel Gaussian normalization',l,d,gaussian_normalization))

    # m_l(z)=(l+1)(1-Phi_d(z))/z^2 <= (l+1)/(2d),
    # from 1-cos x <= x^2/2 and E(theta_1^2)=1/d.
    mmax=arb(l+1)/(2*d)
    alpha=arb(d-2)/2
    # On s=z^2 in [1,2], exact Plancherel substitution is
    # Z = 2/c_d^2 int rho(s)^2/[m_l(sqrt(s))^2 s^alpha] ds.
    # Therefore Z >= K_l ||rho||_2^2 with
    K=arb(2)/(cd*cd*mmax*mmax*(arb(2)**alpha))
    Kclosed=(arb(2)**(arb(d)/2+2))*(arb(d)**2)*((arb(d)/2).gamma()**2)/(arb(l+1)**2)
    if not (K-Kclosed).contains(0):
        raise AssertionError(('physical enstrophy constant simplification',l,K,Kclosed))

    for n in range(1,33):
        cb=math.comb(2*n,n)
        sharp_rho2=arb(2*n+1)*(arb(cb)**2)
        Zfloor=K*sharp_rho2
        # Central-binomial asymptotic comparator: Zfloor/(K*(2/pi)*16^n) -> 1.
        asym=K*(arb(2)/arb.pi())*(arb(16)**n)
        rows.append({
            'l':l,
            'effective_dimension_d':d,
            'hidden_feedback_derivative_orders_0_through':n-1,
            'first_fixed_revealed_dimensionless_derivative_order_n':n,
            'unitary_Hankel_kernel_constant_c_d':str(cd),
            'fixed_band_s_interval':'[1,2]',
            'screened_transfer_upper_bound_m_max':str(mmax),
            'physical_enstrophy_vs_feedback_density_L2_constant_K_l':str(K),
            'sharp_feedback_density_L2_cost_square_for_unit_reveal':str(sharp_rho2),
            'physical_channel_enstrophy_lower_bound_for_unit_reveal':str(Zfloor),
            'enstrophy_floor_divided_by_K_l_16_power_n':str(Zfloor/(K*(arb(16)**n))),
            'ratio_to_asymptotic_K_l_2_over_pi_16_power_n':str(Zfloor/asym),
        })

print(json.dumps({
    'arb_precision_bits':BITS,
    'status':'PASS',
    'cases':len(rows),
    'interpretation':(
        'Normalize the angular toroidal harmonic X_l to unit sphere L2 norm. Then the actual vorticity enstrophy of omega=a(r) r^l X_l is exactly the radial L2 norm int |a|^2 r^(d-1)dr with effective dimension d=2l+3. '
        'Using the self-inverse unitary radial Hankel transform with kernel c_d Phi_d, Plancherel gives Z=int |A(k)|^2 k^(d-1)dk. '
        'For source-scaled spectral variable s=(kL)^2 and feedback density rho, the exact substitution is Z=2/c_d^2 int rho^2/[m_l(sqrt(s))^2 s^((d-2)/2)]ds. '
        'On the fixed physical decay-rate band 1<=s<=2, the spherical-average inequality 1-cos x<=x^2/2 gives m_l<=m_max=(l+1)/(2d). '
        'Combining this exact channel representation with the sharp shifted-Legendre moment theorem yields Z >= K_l |R|^2 (2n+1) C(2n,n)^2 for a fixed nth dimensionless feedback derivative R after hiding all lower derivatives. '
        'Here K_l=2^(d/2+2) d^2 Gamma(d/2)^2/(l+1)^2 in the stated unitary Hankel convention.  Thus for fixed l and fixed band the actual channel enstrophy burden is asymptotically at least K_l (2/pi) 16^n |R|^2. '
        'Restoring a source radius L and vorticity amplitude Omega multiplies the dimensionless channel enstrophy by Omega^2 L^3; the cancellation law itself uses tau=nu t/L^2. '
        'This is a physical fixed-band lower bound. It does not yet exclude shifting or broadening the spectral support, changing source radius, angular mixing, or nonlinear replenishment.'
    ),
    'rows':rows,
},indent=2,allow_nan=False))
