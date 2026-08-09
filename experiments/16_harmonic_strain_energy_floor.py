import json, os, math
import mpmath as mp

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
mp.mp.dps=int(BITS*math.log10(2))+40
pi=mp.pi

# Exact ball-moment formula for even monomials x^(2a)y^(2b)z^(2c) on radius R.
# Integral = 2 R^(2(a+b+c)+3)/(2N+3) * Gamma(a+1/2)Gamma(b+1/2)Gamma(c+1/2)/Gamma(N+3/2)
# (equivalent to radial integral times spherical moment).
def ball_even(a,b,c,R):
    N=a+b+c
    sphere = 2*mp.gamma(a+mp.mpf('0.5'))*mp.gamma(b+mp.mpf('0.5'))*mp.gamma(c+mp.mpf('0.5'))/mp.gamma(N+mp.mpf('1.5'))
    return R**(2*N+3)*sphere/(2*N+3)

# Polynomial represented as dict exponent triple -> coefficient.
def padd(A,B):
    C=dict(A)
    for k,v in B.items(): C[k]=C.get(k,mp.mpf('0'))+v
    return {k:v for k,v in C.items() if v != 0}
def pscale(c,A): return {k:c*v for k,v in A.items()}
def pmul(A,B):
    C={}
    for (a,b,c),u in A.items():
        for (d,e,f),v in B.items():
            k=(a+d,b+e,c+f); C[k]=C.get(k,mp.mpf('0'))+u*v
    return C
def pder(A,j):
    C={}
    for exp,v in A.items():
        e=list(exp)
        if e[j]:
            coeff=v*e[j]; e[j]-=1; C[tuple(e)]=C.get(tuple(e),mp.mpf('0'))+coeff
    return C

def pint_ball(P,R):
    out=mp.mpf('0')
    for (a,b,c),v in P.items():
        if a%2 or b%2 or c%2: continue
        out += v*ball_even(a//2,b//2,c//2,R)
    return out

def grad_energy(phi,R):
    return mp.mpf('0.5')*sum(pint_ball(pmul(pder(phi,j),pder(phi,j)),R) for j in range(3))

# Harmonic solid polynomials: degree 3 and 4 (2D harmonic, hence also 3D harmonic).
H3={(3,0,0):mp.mpf(1),(1,2,0):mp.mpf(-3)}             # x^3-3xy^2
H4={(4,0,0):mp.mpf(1),(2,2,0):mp.mpf(-6),(0,4,0):mp.mpf(1)} # x^4-6x^2y^2+y^4

# Quadratic harmonic potential phi2 = 1/2 x.S.x for symmetric traceless S.
def phi2_from_S(S):
    P={}
    for i in range(3):
        P[tuple(2 if k==i else 0 for k in range(3))]=P.get(tuple(2 if k==i else 0 for k in range(3)),0)+S[i][i]/2
    for i in range(3):
        for j in range(i+1,3):
            e=[0,0,0]; e[i]=1;e[j]=1
            P[tuple(e)]=P.get(tuple(e),0)+S[i][j]
    return P

def frob2(S): return sum(S[i][j]**2 for i in range(3) for j in range(3))

Rs=[mp.mpf(x) for x in ['1e-3','0.1','1','10','1e3']]
matrices=[
    [[mp.mpf('2'),mp.mpf('0'),mp.mpf('0')],[mp.mpf('0'),mp.mpf('-1'),mp.mpf('0')],[mp.mpf('0'),mp.mpf('0'),mp.mpf('-1')]],
    [[mp.mpf('1'),mp.mpf('0.7'),mp.mpf('-0.2')],[mp.mpf('0.7'),mp.mpf('-0.4'),mp.mpf('0.3')],[mp.mpf('-0.2'),mp.mpf('0.3'),mp.mpf('-0.6')]],
    [[mp.mpf('0'),mp.mpf('1'),mp.mpf('0.5')],[mp.mpf('1'),mp.mpf('0'),mp.mpf('-0.25')],[mp.mpf('0.5'),mp.mpf('-0.25'),mp.mpf('0')]],
]
perturbs=[(mp.mpf('0'),mp.mpf('0')),(mp.mpf('0.01'),mp.mpf('0')),(mp.mpf('-0.3'),mp.mpf('0.02')),(mp.mpf('2'),mp.mpf('-0.4'))]
rows=[]
max_floor_rel=mp.mpf('0')
for mi,S in enumerate(matrices):
    tr=sum(S[i][i] for i in range(3))
    if abs(tr)>mp.mpf('1e-70'): raise AssertionError(('S not traceless',mi,tr))
    P2=phi2_from_S(S); sn2=frob2(S)
    for R in Rs:
        floor=(2*pi/15)*sn2*R**5
        E2=grad_energy(P2,R)
        rel=abs(E2-floor)/floor
        max_floor_rel=max(max_floor_rel,rel)
        if rel>mp.mpf('1e-30'): raise AssertionError(('quadratic floor mismatch',mi,R,E2,floor,rel))
        for a,b in perturbs:
            # Scale higher modes by powers of R so coefficients are dimensionless relative perturbations.
            phi=padd(P2,padd(pscale(a/R,H3),pscale(b/(R**2),H4)))
            E=grad_energy(phi,R)
            if E < floor*(1-mp.mpf('1e-28')):
                raise AssertionError(('higher harmonic lowered energy',mi,R,a,b,E,floor))
            rows.append({'matrix':mi,'R':mp.nstr(R,12),'a3':mp.nstr(a,12),'a4':mp.nstr(b,12),
                         'E_over_floor':mp.nstr(E/floor,40)})

# Hodge orthogonality calibration: h=Sx and a smooth toroidal field v=f(r^2) Omega x x.
# v is divergence-free and tangent to every sphere. Its L2 cross term with grad(phi) must vanish.
# We test with f=1+c r^2 and Omega=e_z using exact polynomials.
x={(1,0,0):mp.mpf(1)}; y={(0,1,0):mp.mpf(1)}; z={(0,0,1):mp.mpf(1)}
r2=padd(padd(pmul(x,x),pmul(y,y)),pmul(z,z))
ortho=[]
for c in [mp.mpf('0'),mp.mpf('0.2'),mp.mpf('-0.7')]:
    f=padd({(0,0,0):mp.mpf(1)},pscale(c,r2))
    # Omega=e_z => Omega x x = (-y,x,0)
    v=[pscale(-1,pmul(f,y)),pmul(f,x),{}]
    S=matrices[1]; h=[pder(phi2_from_S(S),j) for j in range(3)]
    crossE=sum(pint_ball(pmul(h[j],v[j]),mp.mpf('1')) for j in range(3))
    if abs(crossE)>mp.mpf('1e-30'): raise AssertionError(('Hodge cross term nonzero',c,crossE))
    Ev=mp.mpf('0.5')*sum(pint_ball(pmul(v[j],v[j]),mp.mpf('1')) for j in range(3))
    Eh=grad_energy(phi2_from_S(S),mp.mpf('1'))
    Et=Eh+Ev+crossE
    if abs(Et-(Eh+Ev))>mp.mpf('1e-30'): raise AssertionError('energy did not split')
    ortho.append({'c':mp.nstr(c,12),'cross_integral':mp.nstr(crossE,20),'Eh':mp.nstr(Eh,30),'Ev':mp.nstr(Ev,30)})

print(json.dumps({
 'precision_bits_requested':BITS,'mpmath_dps':mp.mp.dps,'status':'PASS',
 'exact_floor_constant_2pi_over_15':mp.nstr(2*pi/15,50),
 'max_relative_quadratic_floor_error':mp.nstr(max_floor_rel,20),
 'interpretation':'Inside a ball, the harmonic incompressible component has an exact kinetic-energy floor E_h >= (2pi/15)|S_h(0)|_F^2 d^5. Pure linear strain saturates it; translations/higher harmonic content can only add energy. A divergence-free toroidal component tangent to the sphere is L2-orthogonal to the harmonic component, calibrating the natural potential-vortical energy split.',
 'sample_rows':rows[:12], 'orthogonality_checks':ortho
},indent=2))
