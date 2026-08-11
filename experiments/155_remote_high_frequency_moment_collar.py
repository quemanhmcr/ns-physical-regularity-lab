import json, os, math
from fractions import Fraction as Q
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
one=arb(1); two=arb(2)

def A(q): return arb(q.numerator)/q.denominator

def qmat_inverse(M):
    n=len(M); Aaug=[row[:] + [Q(int(i==j)) for j in range(n)] for i,row in enumerate(M)]
    for c in range(n):
        p=next(i for i in range(c,n) if Aaug[i][c] != 0)
        Aaug[c],Aaug[p]=Aaug[p],Aaug[c]
        piv=Aaug[c][c]; Aaug[c]=[v/piv for v in Aaug[c]]
        for i in range(n):
            if i==c: continue
            f=Aaug[i][c]
            if f: Aaug[i]=[Aaug[i][j]-f*Aaug[c][j] for j in range(2*n)]
    return [row[n:] for row in Aaug]

def moment_int_power(p):
    return Q(2**(p+1)-1,p+1)

def poly_moment(c,q):
    return sum((c[k]*moment_int_power(k+q) for k in range(len(c))),Q(0))

def poly_l2(c):
    return sum((c[i]*c[j]*moment_int_power(i+j) for i in range(len(c)) for j in range(len(c))),Q(0))

def hankel_constants(l):
    d=2*l+3; nu=arb(d)/2-1
    cd=one/((two**nu)*(arb(d)/2).gamma())
    return d,cd

def high_frequency_transfer_lower_constant(l):
    # Phi_d(z)=average cos(z t). For z>=2, use t in [1/(2z),1/z] and its negative.
    d=2*l+3
    sphere_t_density_const=(arb(d)/2).gamma()/(arb.pi().sqrt()*(arb(d-1)/2).gamma())
    c0=one-arb('0.5').cos()
    cstar=c0*sphere_t_density_const*((arb(3)/4)**(arb(d-3)/2))
    # 1-Phi >= cstar/z, hence m_l(z)>=(l+1)cstar/z^3.
    return arb(l+1)*cstar

def low_band_transfer_lower(l):
    d=2*l+3
    # On s=z^2 in [1,2], 1-cos x >= x^2/2-x^4/24 and spherical moments give
    # m_l >= (l+1)(2d+3)/(4d(d+2)).
    return arb(l+1)*arb(2*d+3)/(arb(4*d*(d+2)))

T=arb('1')
rows=[]
for n in range(1,9):
    # Low feedback carrier rho0=1 on [1,2]. Its qth moments are mu_q.
    mu=[moment_int_power(q) for q in range(n)]
    # Dual moment polynomials phi_j on [1,2]: int x^q phi_j dx = delta_qj.
    G=[[moment_int_power(q+k) for k in range(n)] for q in range(n)]
    Gi=qmat_inverse(G)
    for q in range(n):
        for j in range(n):
            chk=sum((G[q][k]*Gi[k][j] for k in range(n)),Q(0))
            if chk != Q(int(q==j)): raise AssertionError(('dual moment inverse',n,q,j,chk))

    M_int=10**(8*n)
    M=Q(M_int)
    # eta_M(M x)=sum_j b_j phi_j(x), b_j=-mu_j M^(-j-1).
    b=[-mu[j]/(M**(j+1)) for j in range(n)]
    coeff=[sum((Gi[k][j]*b[j] for j in range(n)),Q(0)) for k in range(n)]
    for q in range(n):
        corr_moment=(M**(q+1))*poly_moment(coeff,q)
        if corr_moment != -mu[q]: raise AssertionError(('remote collar moment cancellation',n,q,corr_moment,-mu[q]))
    l2_ds=M*poly_l2(coeff)
    if l2_ds<=0: raise AssertionError(('collar L2',n,l2_ds))
    l1_upper=(arb(M_int)*A(l2_ds)).sqrt()
    future_corr_upper=(-arb(M_int)*T).exp()*l1_upper
    F0=((-T).exp()-(-2*T).exp())/T
    future_lower=F0-future_corr_upper
    if not (future_lower.lower()>0): raise AssertionError(('future carrier lost',n,future_lower,future_corr_upper))

    for l in (2,4,8):
        d,cd=hankel_constants(l); alpha=arb(d-2)/2; beta=arb(3)-alpha # =5/2-l
        bhigh=high_frequency_transfer_lower_constant(l)
        if not (bhigh.lower()>0): raise AssertionError(('high frequency transfer constant',l,bhigh))
        # Exact physical spectral enstrophy weight obeys
        # w(s)=2/[c_d^2 m(s)^2 s^alpha] <= 2/(c_d^2 bhigh^2) s^beta on s in [M,2M].
        edge_factor=(arb(2)**beta) if beta>0 else one
        wc=arb(2)/(cd*cd*bhigh*bhigh)
        Zcorr_upper=wc*edge_factor*(arb(M_int)**beta)*A(l2_ds)
        if not (Zcorr_upper.lower()>=0): raise AssertionError(('collar enstrophy upper',n,l,Zcorr_upper))
        blow=low_band_transfer_lower(l)
        Zlow_upper=arb(2)/(cd*cd*blow*blow) # rho0 L2 mass is one and s^alpha>=1.
        Znormalized_upper=(Zlow_upper+Zcorr_upper)/(future_lower*future_lower)
        Zlow_normalized_upper=Zlow_upper/(F0*F0)
        rows.append({
            'l':l,
            'hidden_initial_feedback_moments_count_n':n,
            'remote_collar_band_start_M':str(M_int),
            'remote_collar_band':'[M,2M]',
            'all_total_feedback_moments_0_through_n_minus_1_cancel_exactly':True,
            'collar_unweighted_L2_mass':str(A(l2_ds)),
            'high_frequency_transfer_lower_constant_b_l_in_m_ge_b_s_minus_3_over_2':str(bhigh),
            'physical_collar_enstrophy_upper_bound':str(Zcorr_upper),
            'collar_future_feedback_absolute_upper_at_tau_1':str(future_corr_upper),
            'low_carrier_future_feedback_at_tau_1':str(F0),
            'total_future_feedback_lower_before_normalization':str(future_lower),
            'normalized_total_physical_enstrophy_upper_for_unit_future_feedback':str(Znormalized_upper),
            'low_carrier_only_normalized_enstrophy_upper_reference':str(Zlow_normalized_upper),
            'normalized_upper_over_low_reference':str(Znormalized_upper/Zlow_normalized_upper),
            'collar_scaling_exponent_M_power_3_over_2_minus_l':str(beta),
        })

print(json.dumps({
    'arb_precision_bits':BITS,
    'status':'PASS',
    'cases':len(rows),
    'interpretation':(
        'A fixed-band hiding floor is not global when the viscous spectrum may migrate. Start with a low feedback carrier rho0=1 on s in [1,2]. '
        'For any finite n, exact dual moment polynomials on x in [1,2] build a remote correction collar eta_M on s in [M,2M] whose moments cancel the first n moments of rho0 exactly. '
        'The collar future contribution at tau=1 is exponentially small, bounded by exp(-M) times its L1 mass. '
        'A direct spherical-average lower bound gives 1-Phi_d(z)>=c_d_star/z for z>=2, hence m_l(z)>=b_l z^-3.  Combined with the exact unitary Hankel enstrophy formula, the collar enstrophy is bounded above by C_l,n M^(3/2-l). '
        'For every angular channel l>=2 this exponent is negative, so for each finite n the moment bookkeeping can be pushed to sufficiently remote high frequency with arbitrarily small additional physical enstrophy while leaving the positive-time low-band feedback essentially unchanged. '
        'The explicit runs use M=10^(8n) and l=2,4,8.  This kills any regularity contradiction based only on finitely many initial heat-feedback moments plus total enstrophy when unrestricted spectral migration is allowed. '
        'It does not provide a persistent nonlinear NS escape for free: the remote collar decays on time O(1/M), so maintaining such cancellations over a time interval requires continual nonlinear replenishment.  That maintenance/ancestry question is the next physical bottleneck.'
    ),
    'rows':rows,
},indent=2,allow_nan=False))
