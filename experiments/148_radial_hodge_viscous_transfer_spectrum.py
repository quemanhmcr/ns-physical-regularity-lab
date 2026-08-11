import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
z0=arb(0); one=arb(1); two=arb(2); four=arb(4)

def phi(l,zr):
    # Normalized regular radial Helmholtz mode in effective dimension d=2l+3:
    # Phi_d(x)=0F1(;d/2;-x^2/4), Phi_d(0)=1.
    b=arb(2*l+3)/2
    q=-(zr*zr)/4
    return q.hypgeom_0f1(b)

def direct_screened_C(l,z):
    # Independent evaluation of C_l[Phi_d(z r)] at L=1 by integrating
    # the 0F1 series termwise:
    # int r Phi dr = 1/2 1F2(1;2,b;q)
    # int r^(2l+2) Phi dr = 1/(2l+3) 0F1(;b+1;q).
    b=arb(2*l+3)/2
    q=-(z*z)/4
    i1=q.hypgeom([one],[arb(2),b])/2
    i2=q.hypgeom_0f1(b+1)/arb(2*l+3)
    return -arb(l+1)/arb(2*l+1)*(i1-i2)

def eigen_residual(l,z,r):
    # Differentiate 0F1 analytically and evaluate
    # a''+(2l+2)a'/r+z^2 a.
    b=arb(2*l+3)/2
    q=-(z*z*r*r)/4
    f0=q.hypgeom_0f1(b)
    f1=q.hypgeom_0f1(b+1)
    f2=q.hypgeom_0f1(b+2)
    ap=-(z*z*r)/(2*b)*f1
    app=-(z*z)/(2*b)*f1 + (z**4*r*r)/(4*b*(b+1))*f2
    return app + arb(2*l+2)/r*ap + z*z*f0

z_strings=['1e-8','1e-4','0.1','1','3','10','30','100']
rows=[]
for l in (2,4,6,8,12):
    d=2*l+3
    small_limit=arb(l+1)/(2*d)
    for zs in z_strings:
        z=arb(zs)
        ph=phi(l,z)
        one_minus=one-ph
        if not (one_minus.lower()>0):
            raise AssertionError(('spectral blind spot on sampled real frequency',l,zs,one_minus))
        c_direct=direct_screened_C(l,z)
        c_transfer=-arb(l+1)*one_minus/(z*z)
        err=c_direct-c_transfer
        if not err.contains(0):
            raise AssertionError(('Hodge spectral transfer identity',l,zs,c_direct,c_transfer,err))
        for rs in ('0.17','0.63','1'):
            er=eigen_residual(l,z,arb(rs))
            if not er.contains(0):
                raise AssertionError(('radial Helmholtz eigenmode',l,zs,rs,er))
        m=-c_transfer
        tau=one/(z*z)
        rows.append({
            'l':l,
            'effective_radial_dimension_d':d,
            'dimensionless_frequency_z_kL':zs,
            'Phi_d_z':str(ph),
            'one_minus_Phi_d_z':str(one_minus),
            'screened_transfer_gain_m_l':str(m),
            'viscous_decay_clock_tau_nu_over_L2_over_nu':str(tau),
            'gain_over_decay_clock':str(m/tau),
            'small_frequency_gain_limit':str(small_limit),
            'gain_over_small_frequency_limit':str(m/small_limit),
            'direct_integral_vs_endpoint_transfer_error':str(err),
            'radial_Helmholtz_residual_r017':str(eigen_residual(l,z,arb('0.17'))),
            'radial_Helmholtz_residual_r063':str(eigen_residual(l,z,arb('0.63'))),
            'radial_Helmholtz_residual_r1':str(eigen_residual(l,z,arb('1'))),
        })

# Exact finite-frequency observability control.  If F(t)=sum_j A_j exp(-lambda_j t)
# has its first N time derivatives zero at one time for N distinct lambdas, the
# Vandermonde system forces every A_j=0.  Since m_l(z)>0 analytically for z>0,
# no finite collection of distinct viscous frequencies has a permanent Hodge-blind combination.
vandermonde=[]
for N in range(1,9):
    lam=[(j+1)*(j+1) for j in range(N)]
    det=1
    for i in range(N):
        for j in range(i+1,N):
            det*=lam[j]-lam[i]
    if det==0:
        raise AssertionError(('Vandermonde degeneracy',N))
    vandermonde.append({'number_of_distinct_viscous_frequencies':N,'squared_frequencies_lambda':lam,'Vandermonde_determinant':str(det),'full_rank':True})

print(json.dumps({
    'arb_precision_bits':BITS,
    'status':'PASS',
    'sampled_transfer_cases':len(rows),
    'finite_packet_rank_cases':len(vandermonde),
    'interpretation':(
        'For a toroidal angular Hodge channel, the radial diffusion operator is the radial Laplacian in effective dimension d=2l+3. '
        'Its normalized regular spectral mode Phi_d(z r)=0F1(;d/2;-z^2 r^2/4) satisfies D_l Phi=-z^2 Phi. '
        'Combining that eigenvalue equation with the already validated endpoint identity gives the exact screened transfer '
        'C_l[Phi_d(z r)]=-(l+1)(1-Phi_d(z))/z^2.  The same quantity is independently evaluated from the integrated 0F1 series. '
        'Analytically Phi_d(z) is the spherical average of cos(z theta_1), hence 1-Phi_d(z)>0 for every real z nonzero: there is no radial viscous spectral blind frequency. '
        'For a single generalized frequency the screened gain is m_l(z)=(l+1)(1-Phi_d(z))/z^2 while its dimensionless viscous lifetime is 1/z^2, so high-frequency screening and short viscous lifetime are tied by the same z^2. '
        'The Vandermonde control records the corresponding finite-packet observability: distinct viscous exponentials cannot remain identically Hodge-silent unless every spectral amplitude vanishes.'
    ),
    'rows':rows,
    'finite_packet_observability':vandermonde,
},indent=2,allow_nan=False))
