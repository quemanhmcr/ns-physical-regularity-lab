import json, os, math
from fractions import Fraction as Q
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def poly_mul(a,b):
    c=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c

def poly_pow(a,n):
    r=[Q(1)]
    for _ in range(n): r=poly_mul(r,a)
    return r

def poly_deriv(a,n=1):
    r=a[:]
    for _ in range(n):
        r=[Q(i)*r[i] for i in range(1,len(r))]
    return r

def poly_int_moment(a,q,lo=Q(1),hi=Q(2)):
    s=Q(0)
    for i,c in enumerate(a):
        p=i+q+1
        s+=c*(hi**p-lo**p)/p
    return s

def A(x): return arb(x.numerator)/x.denominator

rows=[]
for N in range(2,17):
    # eta=(s-1)^N (2-s)^N, nonnegative on [1,2], zero to order N at both endpoints.
    eta=poly_mul(poly_pow([Q(-1),Q(1)],N),poly_pow([Q(2),Q(-1)],N))
    rho=poly_deriv(eta,N-1)
    eta_mass=poly_int_moment(eta,0)
    if eta_mass<=0: raise AssertionError(('eta positivity mass',N,eta_mass))
    # Moment cancellations of rho imply zero heat-feedback derivatives at tau=0.
    moments=[poly_int_moment(rho,q) for q in range(N)]
    if any(moments[q]!=0 for q in range(N-1)):
        raise AssertionError(('compact spectral moment cancellation',N,moments))
    expected=(-Q(1) if (N-1)%2 else Q(1))*Q(math.factorial(N-1))*eta_mass
    if moments[N-1]!=expected:
        raise AssertionError(('first moment reveal',N,moments[N-1],expected))
    # rho is compactly supported and square-integrable.  Its exact L2 mass is finite.
    rho2=poly_mul(rho,rho)
    rho_l2=poly_int_moment(rho2,0)
    if rho_l2<=0: raise AssertionError(('rho L2',N,rho_l2))
    # Endpoint traces of rho vanish because eta has zeros of order N and only N-1 derivatives were taken.
    def peval(p,x): return sum((c*x**i for i,c in enumerate(p)),Q(0))
    r1=peval(rho,Q(1)); r2=peval(rho,Q(2))
    if r1!=0 or r2!=0: raise AssertionError(('rho endpoint support regularity',N,r1,r2))
    rows.append({
        'hiding_order_N_minus_1':N-1,
        'eta_polynomial_degree':len(eta)-1,
        'rho_polynomial_degree':len(rho)-1,
        'eta_positive_mass_integral':str(A(eta_mass)),
        'rho_compact_support_L2_mass_integral':str(A(rho_l2)),
        'zero_feedback_derivative_orders':list(range(N-1)),
        'first_nonzero_feedback_derivative_order':N-1,
        'first_nonzero_spectral_moment':str(A(moments[N-1])),
        'first_nonzero_feedback_derivative_positive_magnitude':str(A(Q(math.factorial(N-1))*eta_mass)),
        'rho_endpoint_values':[str(A(r1)),str(A(r2))],
        'feedback_positive_time_factorization':'tau^(N-1) * integral_1^2 exp(-s tau) eta_N(s) ds > 0',
    })

print(json.dumps({
    'arb_precision_bits':BITS,
    'status':'PASS',
    'cases':len(rows),
    'supported_angular_channels':'every l with positive transfer m_l(sqrt(s)) on 1<=s<=2',
    'interpretation':(
        'Finite-energy Hankel spectral states can realize arbitrarily deep finite-jet screened-feedback hiding. '
        'For eta_N=(s-1)^N(2-s)^N on [1,2] and rho_N=d^(N-1)eta_N/ds^(N-1), the first N-1 spectral moments of rho through degree N-2 vanish exactly. '
        'The heat feedback is F_N(tau)=integral exp(-s tau)rho_N(s)ds=tau^(N-1) integral exp(-s tau)eta_N(s)ds, hence it is strictly positive for every tau>0 but has a zero of order N-1 at tau=0. '
        'Because rho_N is compactly supported and L2, and the Hodge transfer multiplier m_l(sqrt(s)) is continuous and strictly positive on [1,2], dividing rho_N by the transfer and the harmless Hankel measure weight produces a compactly supported L2 Hankel spectrum. '
        'By radial Hankel Plancherel this corresponds to a finite-enstrophy radial channel.  Thus arbitrarily deep finite initial hiding is physically compatible with finite enstrophy; what does not survive is exact silence on a positive time interval.'
    ),
    'rows':rows,
},indent=2,allow_nan=False))
