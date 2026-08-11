import json, os, math
from fractions import Fraction as Q
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
one=arb(1); two=arb(2)

def A(q): return arb(q.numerator)/q.denominator

def qinv(M):
    n=len(M); X=[r[:] + [Q(int(i==j)) for j in range(n)] for i,r in enumerate(M)]
    for c in range(n):
        p=next(i for i in range(c,n) if X[i][c])
        X[c],X[p]=X[p],X[c]
        z=X[c][c]; X[c]=[v/z for v in X[c]]
        for i in range(n):
            if i==c: continue
            f=X[i][c]
            if f:X[i]=[X[i][j]-f*X[c][j] for j in range(2*n)]
    return [r[n:] for r in X]

def mpow(p): return Q(2**(p+1)-1,p+1)
def pmoment(c,q): return sum((c[k]*mpow(k+q) for k in range(len(c))),Q(0))
def pl2(c): return sum((c[i]*c[j]*mpow(i+j) for i in range(len(c)) for j in range(len(c))),Q(0))

def collar_l2(n,M):
    mu=[mpow(q) for q in range(n)]
    G=[[mpow(q+k) for k in range(n)] for q in range(n)]
    Gi=qinv(G)
    b=[-mu[j]/(Q(M)**(j+1)) for j in range(n)]
    c=[sum((Gi[k][j]*b[j] for j in range(n)),Q(0)) for k in range(n)]
    for q in range(n):
        if (Q(M)**(q+1))*pmoment(c,q)!=-mu[q]: raise AssertionError(('moment',n,q))
    return Q(M)*pl2(c)

def derivative_L1_bound(l):
    # p_l(t)=C_l(1-t^2)^l on [-1,1].  Bound ||p_l^(l)||_1 by coefficientwise integration.
    C=(arb(l)+arb('0.5')).gamma()/(arb.pi().sqrt()*arb(math.factorial(l)))
    R=Q(0)
    for j in range(l+1):
        if 2*j < l: continue
        m=2*j-l
        coeff=Q(math.comb(l,j)*math.factorial(2*j),math.factorial(m))
        R += abs(coeff)*Q(2,m+1)
    return C*A(R)

def hankel_cd(l):
    d=2*l+3; nu=arb(d)/2-1
    return one/((two**nu)*(arb(d)/2).gamma())

rows=[]
for l in (2,3,4,6,8):
    H=derivative_L1_bound(l)
    zthreshold=(two*H)**(one/l)
    cd=hankel_cd(l)
    for n in range(1,9):
        M=10**(8*n)
        zM=arb(M).sqrt()
        phibound=H/(zM**l)
        if not (phibound.upper() < arb('0.5')):
            raise AssertionError(('high-frequency Phi bound not below half',l,n,zthreshold,zM,phibound))
        # Hence 1-Phi >=1/2 and m_l(sqrt(s)) >= (l+1)/(2s) throughout [M,2M].
        l2=A(collar_l2(n,M))
        gamma=arb('1.5')-l # exponent in enstrophy weight s^(2-alpha)
        # gamma<0 for all l>=2, so max over [M,2M] occurs at M.
        Zupper=arb(8)/(cd*cd*arb(l+1)**2)*(arb(M)**gamma)*l2
        # Palinstrophy has one extra factor s and is <=2M times the enstrophy upper.
        Pupper=arb(2*M)*Zupper
        rows.append({
            'l':l,
            'hidden_moment_count_n':n,
            'remote_scale_M':str(M),
            'spherical_density_lth_derivative_L1_bound_H_l':str(H),
            'frequency_threshold_for_abs_Phi_le_half':str(zthreshold),
            'actual_band_lower_frequency_sqrt_M':str(zM),
            'certified_abs_Phi_upper_at_sqrt_M':str(phibound),
            'uniform_high_frequency_transfer_bound':'m_l(sqrt(s)) >= (l+1)/(2s) on [M,2M]',
            'collar_feedback_density_L2_mass':str(l2),
            'physical_collar_enstrophy_upper_sharp_high_frequency':str(Zupper),
            'physical_collar_palinstrophy_upper':str(Pupper),
            'enstrophy_scaling_exponent_M_power_one_half_minus_l':str(arb('0.5')-l),
            'palinstrophy_scaling_exponent_M_power_three_halves_minus_l':str(arb('1.5')-l),
        })

print(json.dumps({
    'arb_precision_bits':BITS,
    'status':'PASS',
    'cases':len(rows),
    'interpretation':(
        'For effective dimension d=2l+3 the first-coordinate spherical density is proportional to (1-t^2)^l.  Its first l derivatives below order l vanish at both endpoints. '
        'Integrating the spherical cosine average by parts l times gives |Phi_d(z)| <= ||p_l^(l)||_1 z^-l.  Above an explicit finite threshold this is <=1/2, hence the exact Hodge transfer satisfies m_l(z)>=(l+1)/(2 z^2)=(l+1)/(2s). '
        'Applying this sharper high-frequency transfer bound to the exact remote moment collars of module 155 improves the physical enstrophy upper scaling to O(M^(1/2-l)). '
        'The corresponding palinstrophy, which multiplies spectral enstrophy by one extra s and controls pure viscous vorticity dissipation, is O(M^(3/2-l)). '
        'Both exponents are negative for every l>=2.  Thus the same remote collars can have vanishing enstrophy and vanishing viscous palinstrophy/dissipation throughput while cancelling any fixed finite set of feedback moments. '
        'This kills a second static shortcut: viscosity does not necessarily impose a positive maintenance tax on an already-created remote spectral collar. '
        'The unresolved cost is the nonlinear creation and continual routing of those collars, especially under material circulation ancestry and the inward productive cascade.'
    ),
    'rows':rows,
},indent=2,allow_nan=False))
