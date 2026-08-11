import json, os
from fractions import Fraction as Q
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def M(l,k):
    d=2*l+3
    return -Q(l+1,(2*k+2)*(2*k+d))

def capacitary_pair_monomial(l,k):
    # -(l+1)/(d-2) int_0^1 r^(d-1+2k) (r^(2-d)-1) dr
    d=2*l+3
    return -Q(l+1,d-2)*(Q(1,2*k+2)-Q(1,2*k+d))

def Dl_monomial_coeff(l,k):
    # D_l r^(2k) = lambda_k r^(2k-2)
    if k==0:return Q(0)
    d=2*l+3
    return Q(2*k*(2*k+d-2))

def A(q):return arb(q.numerator)/q.denominator

rows=[]
for l in range(1,13):
    d=2*l+3
    # psi=r^(2-d)-1.  Both powers are radial harmonic away from zero.
    p=2-d
    harmonic_coeff=p*(p+d-2)
    if harmonic_coeff!=0:
        raise AssertionError(('capacitary harmonic exponent',l,d,p,harmonic_coeff))
    for k in range(0,13):
        c1=M(l,k); c2=capacitary_pair_monomial(l,k)
        if c1!=c2:
            raise AssertionError(('screen/potential pairing mismatch',l,k,c1,c2))
        if k>=1:
            lam=Dl_monomial_coeff(l,k)
            lhs=lam*M(l,k-1)
            rhs=-Q(l+1)  # (l+1)(a(0)-a(1)) for r^(2k)
            if lhs!=rhs:
                raise AssertionError(('Green endpoint identity monomial',l,k,lhs,rhs))
        rows.append({
            'l':l,
            'effective_dimension_d':d,
            'screen_exponent_d_minus_2':d-2,
            'radial_monomial_power_2k':2*k,
            'Hodge_screen_coefficient':str(A(c1)),
            'capacitary_pairing_coefficient':str(A(c2)),
            'capacitary_potential_radial_harmonic_coefficient':str(harmonic_coeff),
            'Dl_monomial_coefficient':str(A(Dl_monomial_coeff(l,k))),
        })

print(json.dumps({
    'arb_precision_bits':BITS,
    'status':'PASS',
    'cases':len(rows),
    'interpretation':(
        'A three-dimensional toroidal angular channel of degree l carries the radial operator D_l=d2/dr2+(2l+2)r^-1 d/dr, which is the radial Laplacian in effective dimension d=2l+3. '
        'The Hodge screen kernel satisfies r[1-r^(2l+1)]=r^(d-1)[r^(2-d)-1].  Thus C_l is exactly the weighted pairing with the grounded Newtonian capacitary potential psi=r^(2-d)-1, normalized by its radial capacity d-2. '
        'Because psi is D_l-harmonic away from the center, Green identity turns the bulk viscous operator into the center-minus-source-boundary trace C_l[D_l a]=(l+1)(a(0)-a(1)). '
        'The experiment verifies the identity exactly on a spanning monomial family for l=1..12.  The exponent 2l+1, the endpoint law and the z^-2 spectral transfer are therefore different faces of one effective-dimensional potential-theory structure, not independently chosen kernels.'
    ),
    'rows':rows,
},indent=2,allow_nan=False))
