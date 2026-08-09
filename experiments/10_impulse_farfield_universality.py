import json, os, math
import mpmath as mp

BITS = int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be at least 160')
mp.mp.dps = int(BITS*math.log10(2)) + 35
pi = mp.pi
Gamma = mp.mpf('1')

# Closed filament geometries.  The first four have vector area pi e_z.
def circle(t):
    return (mp.cos(t), mp.sin(t), mp.mpf('0'))
def dcircle(t):
    return (-mp.sin(t), mp.cos(t), mp.mpf('0'))

def ellipse(t):
    return (2*mp.cos(t), mp.mpf('0.5')*mp.sin(t), mp.mpf('0'))
def dellipse(t):
    return (-2*mp.sin(t), mp.mpf('0.5')*mp.cos(t), mp.mpf('0'))

radial_a = mp.mpf('0.25')
radial_c = 1/mp.sqrt(1 + radial_a**2/2)  # fixes enclosed area to pi
def radial3(t):
    r = radial_c*(1 + radial_a*mp.cos(3*t))
    return (r*mp.cos(t), r*mp.sin(t), mp.mpf('0'))
def dradial3(t):
    r = radial_c*(1 + radial_a*mp.cos(3*t))
    rp = radial_c*(-3*radial_a*mp.sin(3*t))
    return (rp*mp.cos(t)-r*mp.sin(t), rp*mp.sin(t)+r*mp.cos(t), mp.mpf('0'))

wav_eps = mp.mpf('0.6')
def wavy3d(t):
    return (mp.cos(t), mp.sin(t), wav_eps*mp.sin(2*t))
def dwavy3d(t):
    return (-mp.sin(t), mp.cos(t), 2*wav_eps*mp.cos(2*t))

loops = [
    ('circle', circle, dcircle),
    ('ellipse_area_matched', ellipse, dellipse),
    ('radial_3fold_area_matched', radial3, dradial3),
    ('nonplanar_wavy_area_matched', wavy3d, dwavy3d),
]

def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def add(a,b): return tuple(a[i]+b[i] for i in range(3))
def sub(a,b): return tuple(a[i]-b[i] for i in range(3))
def mul(c,a): return tuple(c*a[i] for i in range(3))
def norm(a): return mp.sqrt(dot(a,a))

def quad_periodic_scalar(f):
    return mp.quad(f, [0, pi/2, pi, 3*pi/2, 2*pi])
def quad_vec(f):
    return tuple(quad_periodic_scalar(lambda t,j=j: f(t)[j]) for j in range(3))

def impulse(X,dX):
    raw = quad_vec(lambda t: cross(X(t), dX(t)))
    return mul(Gamma/2, raw)

def biot_savart(X,dX,x):
    def integrand(t):
        y = X(t); dy=dX(t); r=sub(x,y); rr=norm(r)
        return mul(Gamma/(4*pi*rr**3), cross(dy,r))
    return quad_vec(integrand)

def dipole_velocity(I,x):
    r=norm(x); e=mul(1/r,x); Ie=dot(I,e)
    return mul(1/(4*pi*r**3), sub(mul(3*Ie,e), I))

# Precision-scaled comparison tolerance: quadrature is arbitrary precision but not interval-enclosed.
tol_imp = mp.mpf(10) ** (-min(35, mp.mp.dps//3))
Ds = [mp.mpf('20'), mp.mpf('50'), mp.mpf('100')]
directions = [
    ('axial', (mp.mpf('0'),mp.mpf('0'),mp.mpf('1'))),
    ('equatorial', (mp.mpf('1'),mp.mpf('0'),mp.mpf('0'))),
    ('oblique', mul(1/mp.sqrt(3),(mp.mpf('1'),mp.mpf('1'),mp.mpf('1')))),
]
rows=[]
last_errors=[]
for name,X,dX in loops:
    I=impulse(X,dX)
    targetI=(mp.mpf('0'),mp.mpf('0'),pi)
    imp_rel=norm(sub(I,targetI))/norm(targetI)
    if imp_rel > tol_imp:
        raise AssertionError((name,'vector-area/impulse mismatch',imp_rel,I))
    for dname,e in directions:
        errs=[]
        for D in Ds:
            x=mul(D,e)
            u=biot_savart(X,dX,x)
            ud=dipole_velocity(I,x)
            rel=norm(sub(u,ud))/norm(ud)
            errs.append(rel)
        # Far-field asymptotics must improve as observation distance grows.
        if not (errs[2] < errs[1] and errs[1] < errs[0]):
            raise AssertionError((name,dname,'far-field error not monotonically decreasing',errs))
        if errs[-1] > mp.mpf('0.02'):
            raise AssertionError((name,dname,'dipole universality not reached by D=100',errs[-1]))
        last_errors.append(errs[-1])
        rows.append({
            'loop': name,
            'direction': dname,
            'impulse': [mp.nstr(v,30) for v in I],
            'relative_errors_D20_D50_D100': [mp.nstr(v,24) for v in errs],
        })

# Universality gate: all microscopic shapes with the same impulse converge to the same far field.
max_D100=max(last_errors)
print(json.dumps({
    'precision_bits_requested':BITS,
    'mpmath_dps':mp.mp.dps,
    'status':'PASS',
    'loops':len(loops),
    'directions':len(directions),
    'max_relative_error_at_D100':mp.nstr(max_D100,30),
    'interpretation':'Closed vortex geometry cancels the monopole exactly. Distinct planar/nonplanar loops with the same hydrodynamic impulse converge to the same 1/r^3 dipole velocity field; microscopic shape is forgotten in the far field.',
    'rows':rows,
},indent=2))
