import json, os, math
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
one=arb(1)

def phi(l,z):
    b=arb(2*l+3)/2
    return (-(z*z)/4).hypgeom_0f1(b)

def gain(l,z):
    return arb(l+1)*(one-phi(l,z))/(z*z)

rows=[]
for l in (2,4,8,12):
    for N in range(2,13):
        lam=[j+1 for j in range(N)]
        A=[((-1)**j)*math.comb(N-1,j) for j in range(N)]
        # Exact finite-difference moment cancellations.
        moments=[]
        for q in range(N):
            s=sum(A[j]*(lam[j]**q) for j in range(N))
            moments.append(s)
        if any(moments[q]!=0 for q in range(N-1)):
            raise AssertionError(('finite spectral jet cancellation',l,N,moments))
        expected=(((-1)**(N-1))*math.factorial(N-1))
        if moments[N-1]!=expected:
            raise AssertionError(('first revealed derivative',l,N,moments[N-1],expected))

        gains=[]; coeffs=[]
        for j in range(N):
            z=arb(lam[j]).sqrt()
            m=gain(l,z)
            if not (m.lower()>0):
                raise AssertionError(('nonpositive individual Hodge gain',l,N,j,m))
            gains.append(m)
            # C_j=-m_j c_j; choose c_j=-A_j/m_j so total screened signal is sum A_j e^-lambda_j tau.
            coeffs.append(-arb(A[j])/m)

        sample_rows=[]
        for ts in ('1e-6','0.001','0.01','0.1','1'):
            t=arb(ts)
            direct=sum((arb(A[j])*(-arb(lam[j])*t).exp() for j in range(N)),arb(0))
            x=(-t).exp()
            closed=x*(one-x)**(N-1)
            err=direct-closed
            if not err.contains(0):
                raise AssertionError(('spectral packet closed form',l,N,ts,direct,closed,err))
            if not (closed.lower()>0):
                raise AssertionError(('positive-time packet unexpectedly silent',l,N,ts,closed))
            sample_rows.append({'tau':ts,'feedback_signal':str(closed),'direct_vs_closed_error':str(err)})

        tpeak=arb(N).log()
        xpeak=(-tpeak).exp()
        peak=xpeak*(one-xpeak)**(N-1)
        exact_peak=arb(1)/N*(arb(N-1)/N)**(N-1)
        if not (peak-exact_peak).contains(0):
            raise AssertionError(('peak identity',l,N,peak,exact_peak))
        rows.append({
            'l':l,
            'number_of_visible_viscous_frequencies_N':N,
            'lambda_j':lam,
            'screened_packet_weights_A_j':A,
            'individual_screened_gains_m_j':[str(v) for v in gains],
            'generalized_spectral_vorticity_coefficients_c_j':[str(v) for v in coeffs],
            'exact_zero_feedback_derivative_orders':list(range(N-1)),
            'first_nonzero_feedback_derivative_order':N-1,
            'first_nonzero_feedback_derivative':str(expected),
            'positive_time_closed_form':'exp(-tau)*(1-exp(-tau))^(N-1)',
            'peak_dimensionless_time_log_N':str(tpeak),
            'peak_feedback_signal':str(peak),
            'time_samples':sample_rows,
        })

print(json.dumps({
    'arb_precision_bits':BITS,
    'status':'PASS',
    'cases':len(rows),
    'interpretation':(
        'The absence of individual radial viscous blind frequencies does not forbid signed interference between visible frequencies. '
        'For lambda_j=j+1 and screened packet weights A_j=(-1)^j binomial(N-1,j), the exact heat feedback is F_N(tau)=exp(-tau)(1-exp(-tau))^(N-1). '
        'Its first N-1 Taylor coefficients through derivative order N-2 vanish exactly, while F_N(tau)>0 for every positive tau and peaks at tau=log N. '
        'Thus arbitrarily deep finite-jet hiding survives in the true viscous spectral language; what is forbidden is exact silence on a positive time interval for a finite packet. '
        'The reported c_j are generalized spectral coefficients needed to compensate the strictly positive individual Hodge gains.  They are not interpreted as a finite-energy norm because discrete radial Helmholtz modes are generalized, not L2, eigenfunctions.'
    ),
    'rows':rows,
},indent=2,allow_nan=False))
