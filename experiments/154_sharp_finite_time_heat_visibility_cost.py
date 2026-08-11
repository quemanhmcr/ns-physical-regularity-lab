import json, os, math
from fractions import Fraction as Q
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
one=arb(1); two=arb(2)

def shifted_legendre_coeffs(n):
    return [Q(((-1)**(n-k))*math.comb(n,k)*math.comb(n+k,k)) for k in range(n+1)]

def A(q): return arb(q.numerator)/q.denominator

def monomial_heat_integral(k,T):
    # int_0^1 x^k exp(-T x) dx = 1/(k+1) 1F1(k+1;k+2;-T)
    return (-T).hypgeom([arb(k+1)],[arb(k+2)])/arb(k+1)

def legendre_heat_integral_poly(j,T):
    c=shifted_legendre_coeffs(j)
    return sum((A(c[k])*monomial_heat_integral(k,T) for k in range(j+1)),arb(0))

def sph_i(j,a):
    # modified spherical Bessel i_j(a)=sqrt(pi/(2a)) I_{j+1/2}(a)
    return (arb.pi()/(2*a)).sqrt()*a.bessel_i(arb(j)+arb('0.5'))

def legendre_heat_integral_bessel(j,T):
    a=T/2
    v=(-T/2).exp()*sph_i(j,a)
    return -v if (j%2) else v

def positive_tail_bounds(n,T,extra=28):
    # ||(I-P_{n-1}) exp(-T x)||_L2(0,1)^2
    # = sum_{j=n}^inf (2j+1) I_j(T)^2, all terms positive.
    M=n+extra
    partial=arb(0)
    for j in range(n,M+1):
        I=legendre_heat_integral_bessel(j,T)
        partial += arb(2*j+1)*I*I
    J=M+1
    a=T/2
    # Series upper bound for i_J(a). First hypergeometric/Bessel series term times geometric correction.
    B=(arb.pi().sqrt()*(a**J))/(arb(2)**(J+1)*(arb(J)+arb('1.5')).gamma())
    r=(a*a)/(4*(arb(J)+arb('1.5')))
    if not (r.upper()<1): raise AssertionError(('Bessel series ratio',n,T,J,r))
    U=B/(1-r)
    qratio=(arb(2*J+3)/arb(2*J+1))*(a/(2*(arb(J)+arb('1.5'))))**2
    if not (qratio.upper()<1): raise AssertionError(('j-tail ratio',n,T,J,qratio))
    rem=arb(2*J+1)*U*U/(1-qratio)
    lower=arb(partial.lower())
    upper=arb((partial+rem).upper())
    if not (lower>0 and upper>=lower): raise AssertionError(('positive tail bracket',n,T,lower,upper))
    return lower,upper,rem

def physical_K(l):
    d=2*l+3
    nu=arb(d)/2-1
    cd=one/((two**nu)*(arb(d)/2).gamma())
    mmax=arb(l+1)/(2*d)
    alpha=arb(d-2)/2
    return arb(2)/(cd*cd*mmax*mmax*(two**alpha))

rows=[]
Tvals=['0.05','0.1','0.25','0.5','1','2','4']
for Ts in Tvals:
    T=arb(Ts)
    for n in range(1,17):
        # Independent polynomial-vs-Bessel calibration for the first several Legendre channels.
        calerr=arb(0)
        for j in range(min(n+2,8)):
            e=legendre_heat_integral_poly(j,T)-legendre_heat_integral_bessel(j,T)
            calerr += e*e
        if not calerr.contains(0): raise AssertionError(('Legendre heat observer calibration',Ts,n,calerr))
        qlo,qhi,rem=positive_tail_bounds(n,T)
        # On the physical band s in [1,2], exp(-sT)=exp(-T) exp(-Tx).
        qlo_s=(-2*T).exp()*qlo
        qhi_s=(-2*T).exp()*qhi
        # Minimal ||rho|| for F(T)=1 and rho orthogonal to polynomials degree<n.
        cost_lo=one/qhi_s.sqrt()  # rigorous lower bound on the exact sharp cost
        cost_hi=one/qlo_s.sqrt()
        if not (cost_hi>=cost_lo): raise AssertionError(('cost bracket',Ts,n,cost_lo,cost_hi))
        # small-T leading Riesz asymptotic from the first surviving Legendre component
        asym=(T**n)/(arb(math.factorial(n))*arb(math.comb(2*n,n))*arb(2*n+1).sqrt())
        qnorm_lo=qlo_s.sqrt(); qnorm_hi=qhi_s.sqrt()
        rows.append({
            'target_dimensionless_time_tau_star':Ts,
            'hidden_initial_feedback_derivative_orders_0_through':n-1,
            'Legendre_tail_norm_square_lower':str(qlo_s),
            'Legendre_tail_norm_square_upper':str(qhi_s),
            'positive_tail_remainder_bound_before_exp_shift':str(rem),
            'sharp_feedback_density_L2_cost_lower_for_unit_F_at_tau':str(cost_lo),
            'sharp_feedback_density_L2_cost_upper_for_unit_F_at_tau':str(cost_hi),
            'cost_bracket_relative_width':str((cost_hi-cost_lo)/cost_lo),
            'tail_norm_over_small_time_leading_scale_lower':str(qnorm_lo/asym),
            'tail_norm_over_small_time_leading_scale_upper':str(qnorm_hi/asym),
            'observer_calibration_error_square':str(calerr),
            'physical_enstrophy_floor_l2':str(physical_K(2)*cost_lo*cost_lo),
            'physical_enstrophy_floor_l4':str(physical_K(4)*cost_lo*cost_lo),
            'physical_enstrophy_floor_l8':str(physical_K(8)*cost_lo*cost_lo),
        })

print(json.dumps({
    'arb_precision_bits':BITS,
    'status':'PASS',
    'cases':len(rows),
    'interpretation':(
        'On the fixed viscous band s in [1,2], impose exact initial feedback-moment cancellations through degree n-1 and require the screened heat feedback to equal one at a specified positive time tau_star. '
        'In L2(1,2), the Riesz-minimal feedback density is proportional to the component of exp(-s tau_star) orthogonal to polynomials of degree below n. '
        'Its exact sharp cost is the reciprocal of that Legendre-tail norm.  The experiment evaluates the tail as a positive modified-spherical-Bessel series, with a rigorous analytic remainder bound, so deep hiding is not inferred by subtracting nearly equal parent norms. '
        'Independent polynomial-hypergeometric and Bessel formulas for each Legendre heat coefficient are cross-checked. '
        'Multiplying the rigorous sharp-cost lower bound by the validated physical Hankel constant K_l gives an actual channel-enstrophy lower bound for producing unit screened feedback at tau_star after n hidden initial derivatives. '
        'This is a fixed-band finite-time observability cost.  Moving or widening the spectral band and nonlinear replenishment remain open escapes.'
    ),
    'rows':rows,
},indent=2,allow_nan=False))
