import json, os, math
import mpmath as mp

BITS = int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be at least 160')
mp.mp.dps = int(BITS*math.log10(2)) + 35
pi=mp.pi

# ---------- vector utilities ----------
def add(a,b): return tuple(a[i]+b[i] for i in range(3))
def sub(a,b): return tuple(a[i]-b[i] for i in range(3))
def mul(c,a): return tuple(c*a[i] for i in range(3))
def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def norm(a): return mp.sqrt(dot(a,a))
def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def quad_periodic_scalar(f):
    return mp.quad(f,[0,pi/2,pi,3*pi/2,2*pi])
def quad_vec(f):
    return tuple(quad_periodic_scalar(lambda t,j=j:f(t)[j]) for j in range(3))

# ---------- exact circular-ring ladder on axis ----------
R=mp.mpf('1'); Gamma=mp.mpf('1'); a=mp.mpf('3')

def ring_axis(z,z0,sgn):
    return sgn*Gamma*R**2/(2*(R**2+(z-z0)**2)**mp.mpf('1.5'))

def stack_velocity(z,stack):
    return sum(ring_axis(z,z0,s) for z0,s in stack)

single=[(mp.mpf('0'),mp.mpf('1'))]
pair=[(a/2,mp.mpf('1')),(-a/2,mp.mpf('-1'))]
quartet=[(-3*a/2,mp.mpf('1')),(-a/2,mp.mpf('-1')),(a/2,mp.mpf('-1')),(3*a/2,mp.mpf('1'))]

# Moment checks for the signed impulse carriers (common impulse magnitude pi Gamma R^2).
def signed_moments(stack,max_order=3):
    return [sum(s*z0**n for z0,s in stack) for n in range(max_order+1)]

m_single=signed_moments(single)
m_pair=signed_moments(pair)
m_quartet=signed_moments(quartet)
if abs(m_pair[0]) != 0: raise AssertionError(('pair net impulse sign sum not zero',m_pair))
if abs(m_quartet[0]) != 0 or abs(m_quartet[1]) != 0:
    raise AssertionError(('quartet did not cancel first two moment levels',m_quartet))

expected_constants={
    'single': Gamma*R**2/2,                       # u z^3
    'pair': 3*Gamma*R**2*a/2,                   # u z^4
    'quartet': 12*Gamma*R**2*a**2,              # u z^5
}
stacks=[('single',single,3),('pair',pair,4),('quartet',quartet,5)]
Zs=[mp.mpf('1e2'),mp.mpf('1e4'),mp.mpf('1e6')]
ladder=[]
for name,stack,power in stacks:
    vals=[]; errs=[]
    C=expected_constants[name]
    for z in Zs:
        u=stack_velocity(z,stack)
        scaled=u*z**power
        rel=abs(scaled-C)/abs(C)
        vals.append((z,u,scaled)); errs.append(rel)
    if not (errs[2] < errs[1] and errs[1] < errs[0]):
        raise AssertionError((name,'scaled asymptotic error not decreasing',errs))
    if errs[-1] > mp.mpf('1e-9'):
        raise AssertionError((name,'asymptotic constant not reached',errs[-1]))
    # Local observed power between last two radii.
    u1=abs(vals[-2][1]); u2=abs(vals[-1][1])
    observed=-mp.log(u2/u1)/mp.log(Zs[-1]/Zs[-2])
    if abs(observed-power) > mp.mpf('1e-6'):
        raise AssertionError((name,'wrong far-field power',observed,power))
    ladder.append({
        'name':name,
        'predicted_velocity_power':power,
        'observed_power_last_pair':mp.nstr(observed,40),
        'signed_source_moments_n0_n1_n2_n3':[mp.nstr(v,30) for v in signed_moments(stack,3)],
        'expected_scaled_constant':mp.nstr(C,40),
        'scaled_values_z1e2_z1e4_z1e6':[mp.nstr(v[2],40) for v in vals],
        'relative_errors':[mp.nstr(e,24) for e in errs],
    })

# ---------- microscopic-shape universality for the impulse-neutral pair ----------
def circle(t): return (mp.cos(t),mp.sin(t),mp.mpf('0'))
def dcircle(t): return (-mp.sin(t),mp.cos(t),mp.mpf('0'))
def ellipse(t): return (2*mp.cos(t),mp.mpf('0.5')*mp.sin(t),mp.mpf('0'))
def dellipse(t): return (-2*mp.sin(t),mp.mpf('0.5')*mp.cos(t),mp.mpf('0'))
radial_a=mp.mpf('0.25'); radial_c=1/mp.sqrt(1+radial_a**2/2)
def radial3(t):
    rr=radial_c*(1+radial_a*mp.cos(3*t)); return (rr*mp.cos(t),rr*mp.sin(t),mp.mpf('0'))
def dradial3(t):
    rr=radial_c*(1+radial_a*mp.cos(3*t)); rp=-3*radial_c*radial_a*mp.sin(3*t)
    return (rp*mp.cos(t)-rr*mp.sin(t),rp*mp.sin(t)+rr*mp.cos(t),mp.mpf('0'))
wav=mp.mpf('0.6')
def wavy(t): return (mp.cos(t),mp.sin(t),wav*mp.sin(2*t))
def dwavy(t): return (-mp.sin(t),mp.cos(t),2*wav*mp.cos(2*t))
loops=[('circle',circle,dcircle),('ellipse',ellipse,dellipse),('radial3',radial3,dradial3),('wavy3d',wavy,dwavy)]

def impulse(X,dX,sgn=1):
    raw=quad_vec(lambda t:cross(X(t),dX(t)))
    return mul(mp.mpf(sgn)/2,raw)

def bs_shifted(X,dX,center,sgn,x):
    def integrand(t):
        y=add(X(t),center); dy=dX(t); r=sub(x,y); rr=norm(r)
        return mul(mp.mpf(sgn)/(4*pi*rr**3),cross(dy,r))
    return quad_vec(integrand)

def dipole_velocity(I,x):
    r=norm(x); e=mul(1/r,x); Ie=dot(I,e)
    return mul(1/(4*pi*r**3),sub(mul(3*Ie,e),I))

def neutral_pair_leading(I,avec,x):
    # - (a.grad) of the impulse dipole field.
    r=norm(x); e=mul(1/r,x)
    Ie=dot(I,e); ae=dot(avec,e); Ia=dot(I,avec)
    bracket=sub(sub(sub(mul(5*Ie*ae,e),mul(Ie,avec)),mul(Ia,e)),mul(ae,I))
    return mul(3/(4*pi*r**4),bracket)

sep=mp.mpf('3'); cplus=(sep/2,mp.mpf('0'),mp.mpf('0')); cminus=(-sep/2,mp.mpf('0'),mp.mpf('0'))
avec=sub(cplus,cminus)
targetI=(mp.mpf('0'),mp.mpf('0'),pi)
directions=[
    ('generic1',mul(1/mp.sqrt(3),(mp.mpf('1'),mp.mpf('1'),mp.mpf('1')))),
    ('generic2',mul(1/mp.sqrt(14),(mp.mpf('2'),mp.mpf('-1'),mp.mpf('3')))),
]
Ds=[mp.mpf('50'),mp.mpf('100'),mp.mpf('200')]
tol_imp=mp.mpf(10)**(-min(35,mp.mp.dps//3))
shape_rows=[]; final_errors=[]
for name,X,dX in loops:
    I=impulse(X,dX,1)
    if norm(sub(I,targetI))/norm(targetI)>tol_imp:
        raise AssertionError((name,'shape impulse mismatch',I))
    for dname,e in directions:
        errs=[]
        for D in Ds:
            x=mul(D,e)
            u=add(bs_shifted(X,dX,cplus,1,x),bs_shifted(X,dX,cminus,-1,x))
            uq=neutral_pair_leading(I,avec,x)
            rel=norm(sub(u,uq))/norm(uq)
            errs.append(rel)
        if not (errs[2]<errs[1] and errs[1]<errs[0]):
            raise AssertionError((name,dname,'neutral-pair far error not decreasing',errs))
        if errs[-1]>mp.mpf('0.06'):
            raise AssertionError((name,dname,'shape universality not reached by D=200',errs[-1]))
        final_errors.append(errs[-1])
        shape_rows.append({
            'loop':name,'direction':dname,
            'relative_errors_D50_D100_D200':[mp.nstr(v,24) for v in errs],
        })

print(json.dumps({
    'precision_bits_requested':BITS,
    'mpmath_dps':mp.mp.dps,
    'status':'PASS',
    'ladder_interpretation':'Exact closed circular donor stacks show velocity decay r^-3, r^-4, r^-5 as successive physically generated impulse moments are cancelled.',
    'pair_universality_interpretation':'Distinct closed-loop microscopic geometries with the same opposite impulses and separation converge to the same r^-4 impulse-separation field.',
    'max_pair_shape_relative_error_at_D200':mp.nstr(max(final_errors),30),
    'ladder':ladder,
    'shape_rows':shape_rows,
},indent=2))
