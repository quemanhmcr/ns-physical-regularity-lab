import json, os, math
from fractions import Fraction as Q
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def shifted_legendre_coeffs(n):
    # L_n(x)=P_n(2x-1)=sum_k (-1)^(n-k) C(n,k) C(n+k,k) x^k.
    return [Q(((-1)**(n-k))*math.comb(n,k)*math.comb(n+k,k)) for k in range(n+1)]

def poly_mul(a,b):
    c=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c

def int01_poly(a):
    return sum((c/Q(i+1) for i,c in enumerate(a)),Q(0))

def moment_s(a,q):
    # integral_0^1 (1+x)^q a(x) dx, corresponding to s=x+1 in [1,2].
    total=Q(0)
    for j in range(q+1):
        b=Q(math.comb(q,j))
        total += b*sum((c/Q(i+j+1) for i,c in enumerate(a)),Q(0))
    return total

def eta_coeffs(N):
    # x^N(1-x)^N
    a=[Q(0)]*N+[Q(1)]
    b=[Q(((-1)**k)*math.comb(N,k)) for k in range(N+1)]
    return poly_mul(a,b)

def deriv(a,n):
    r=a[:]
    for _ in range(n): r=[Q(i)*r[i] for i in range(1,len(r))]
    return r

def A(q): return arb(q.numerator)/q.denominator

rows=[]
sharp_norms=[]
for n in range(1,33):
    L=shifted_legendre_coeffs(n)
    lead=L[-1]
    expected_lead=Q(math.comb(2*n,n))
    if lead!=expected_lead:
        raise AssertionError(('Legendre leading coefficient',n,lead,expected_lead))
    moms=[moment_s(L,q) for q in range(n+1)]
    if any(moms[q]!=0 for q in range(n)):
        raise AssertionError(('Legendre moment orthogonality',n,moms))
    expected_mn=Q(1,(2*n+1)*math.comb(2*n,n))
    if moms[n]!=expected_mn:
        raise AssertionError(('Legendre first visible moment',n,moms[n],expected_mn))
    norm2=int01_poly(poly_mul(L,L))
    if norm2!=Q(1,2*n+1):
        raise AssertionError(('Legendre L2 norm',n,norm2))
    # Normalize nth moment to +1. This is the Riesz-minimal L2 density.
    scale=Q(1)/moms[n]
    sharp_norm2=scale*scale*norm2
    expected_sharp=Q(2*n+1)*Q(math.comb(2*n,n)**2)
    if sharp_norm2!=expected_sharp:
        raise AssertionError(('sharp cost closed form',n,sharp_norm2,expected_sharp))
    sharp_norm=A(sharp_norm2).sqrt()

    # Compare the compact-support smooth-ish bump construction from module 151 at N=n+1.
    N=n+1
    eta=eta_coeffs(N)
    rho=deriv(eta,n)
    rho2=int01_poly(poly_mul(rho,rho))
    r_mn=moment_s(rho,n)
    if r_mn==0: raise AssertionError(('bump reveal moment',n))
    bump_cost2=rho2/(r_mn*r_mn)
    if bump_cost2 < sharp_norm2:
        raise AssertionError(('Riesz lower bound violated',n,bump_cost2,sharp_norm2))
    bump_over_sharp=(A(bump_cost2)/A(sharp_norm2)).sqrt()
    sharp_norms.append(sharp_norm)
    rows.append({
        'hidden_derivative_orders_0_through':n-1,
        'first_fixed_revealed_derivative_order_n':n,
        'central_binomial_C_2n_n':str(math.comb(2*n,n)),
        'shifted_Legendre_L2_norm_square':str(A(norm2)),
        'shifted_Legendre_nth_moment':str(A(moms[n])),
        'sharp_minimum_L2_cost_square_for_unit_nth_moment':str(A(sharp_norm2)),
        'sharp_minimum_L2_cost_for_unit_nth_moment':str(sharp_norm),
        'successive_cost_ratio_to_previous':None,
        'cost_divided_by_4_power_n':str(sharp_norm/(arb(4)**n)),
        'compact_bump_cost_over_sharp_minimum':str(bump_over_sharp),
    })
for i in range(1,len(rows)):
    rows[i]['successive_cost_ratio_to_previous']=str(sharp_norms[i]/sharp_norms[i-1])

print(json.dumps({
    'arb_precision_bits':BITS,
    'status':'PASS',
    'cases':len(rows),
    'interpretation':(
        'On the fixed viscous decay-rate band s in [1,2], hiding feedback derivatives 0 through n-1 means the feedback spectral density rho is orthogonal in L2 to every polynomial of degree below n. '
        'For a fixed nth spectral moment R, the Riesz-minimal density is the shifted Legendre polynomial of degree n. '
        'Its leading coefficient is C(2n,n), its L2 norm square is 1/(2n+1), and its nth moment is 1/[(2n+1)C(2n,n)]. '
        'Therefore the sharp minimum L2 burden for unit nth revealed moment is sqrt(2n+1) C(2n,n), asymptotic to sqrt(2/pi) 4^n. '
        'The compact bump construction of module 151 is checked to lie above this sharp minimum. '
        'This is an unavoidable fixed-band cancellation cost, not yet a global Navier-Stokes cost: moving or broadening the spectral support remains an escape to attack separately.'
    ),
    'rows':rows,
},indent=2,allow_nan=False))
