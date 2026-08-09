import json, os, math
import mpmath as mp

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
mp.mp.dps=int(BITS*math.log10(2))+40
pi=mp.pi

# Sphere average monomials. Integral over unit sphere n_x^(2a)n_y^(2b)n_z^(2c) dOmega.
def sphere_even(a,b,c):
    return 2*mp.gamma(a+mp.mpf('0.5'))*mp.gamma(b+mp.mpf('0.5'))*mp.gamma(c+mp.mpf('0.5'))/mp.gamma(a+b+c+mp.mpf('1.5'))

def sphere_poly_integral(P):
    out=mp.mpf('0')
    for (a,b,c),v in P.items():
        if a%2 or b%2 or c%2: continue
        out += v*sphere_even(a//2,b//2,c//2)
    return out

def padd(A,B):
    C=dict(A)
    for k,v in B.items(): C[k]=C.get(k,mp.mpf('0'))+v
    return {k:v for k,v in C.items() if v!=0}
def pscale(c,A): return {k:c*v for k,v in A.items()}
def pmul(A,B):
    C={}
    for e,u in A.items():
        for f,v in B.items():
            k=tuple(e[i]+f[i] for i in range(3)); C[k]=C.get(k,mp.mpf('0'))+u*v
    return C

nx={(1,0,0):mp.mpf(1)}; ny={(0,1,0):mp.mpf(1)}; nz={(0,0,1):mp.mpf(1)}
ns=[nx,ny,nz]

Ss=[
 [[mp.mpf('2'),0,0],[0,mp.mpf('-1'),0],[0,0,mp.mpf('-1')]],
 [[mp.mpf('1'),mp.mpf('0.7'),mp.mpf('-0.2')],[mp.mpf('0.7'),mp.mpf('-0.4'),mp.mpf('0.3')],[mp.mpf('-0.2'),mp.mpf('0.3'),mp.mpf('-0.6')]],
]
Ss=[[[mp.mpf(v) for v in row] for row in S] for S in Ss]

def q_on_sphere(S):
    q={}
    for i in range(3):
        for j in range(3):
            q=padd(q,pscale(S[i][j],pmul(ns[i],ns[j])))
    return q

def recover_S_from_g(g,R):
    T=[[mp.mpf('0') for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            nij=pmul(ns[i],ns[j])
            if i==j: nij=padd(nij,{(0,0,0):mp.mpf('-1')/3})
            T[i][j]=(15/(8*pi*R))*sphere_poly_integral(pmul(g,nij))
    return T

def ferr(A,B):
    num=mp.sqrt(sum((A[i][j]-B[i][j])**2 for i in range(3) for j in range(3)))
    den=max(mp.sqrt(sum(B[i][j]**2 for i in range(3) for j in range(3))),mp.mpf('1e-90'))
    return num/den

# Higher normal-flux harmonics that should be invisible to S_h at the center.
# degree 1 translation: n_x; degree 3: n_x^3-3n_x n_y^2; degree 4 harmonic in x,y.
G1=nx
G3=padd(pmul(pmul(nx,nx),nx),pscale(-3,pmul(nx,pmul(ny,ny))))
G4=padd(padd(pmul(pmul(nx,nx),pmul(nx,nx)),pscale(-6,pmul(pmul(nx,nx),pmul(ny,ny)))),pmul(pmul(ny,ny),pmul(ny,ny)))

projection_rows=[]
for si,S in enumerate(Ss):
    q=q_on_sphere(S)
    for R in [mp.mpf(x) for x in ['1e-4','0.1','1','10','1e4']]:
        # pure quadratic normal flux g=R n.S.n, plus unrelated degrees 1,3,4.
        for a1,a3,a4 in [(0,0,0),('2','0.3','-0.2'),('-5','2','3')]:
            g=pscale(R,q)
            g=padd(g,pscale(mp.mpf(a1),G1)); g=padd(g,pscale(mp.mpf(a3),G3)); g=padd(g,pscale(mp.mpf(a4),G4))
            Srec=recover_S_from_g(g,R)
            er=ferr(Srec,S)
            if er>mp.mpf('1e-28'): raise AssertionError(('boundary projector contaminated',si,R,a1,a3,a4,er,Srec,S))
            projection_rows.append({'S':si,'R':mp.nstr(R,12),'extras':[a1,a3,a4],'relative_error':mp.nstr(er,20)})

# Exact self-contained vortical carrier in unit ball:
# v=(1-5r^2/3) Sx + (2/3)(x.S.x)x.
# It is divergence-free and v.n=0 at r=1. On a sub-sphere radius R,
# v.n = R(1-R^2)(n.S.n), so the Hodge harmonic center strain is (1-R^2)S.
carrier=[]
for si,S in enumerate(Ss):
    q=q_on_sphere(S)
    for R in [mp.mpf(x) for x in ['0.001','0.01','0.1','0.3','0.5','0.7071067811865475244','0.9','1.0']]:
        g=pscale(R*(1-R**2),q)
        Srec=recover_S_from_g(g,R)
        Starget=[[ (1-R**2)*S[i][j] for j in range(3)] for i in range(3)]
        er=ferr(Srec,Starget) if R<1 else mp.sqrt(sum(Srec[i][j]**2 for i in range(3) for j in range(3)))
        if er>mp.mpf('1e-27'): raise AssertionError(('carrier Hodge profile mismatch',si,R,er))
        harmonic_fraction=1-R**2
        vortical_fraction=R**2
        carrier.append({'S':si,'R_over_carrier':mp.nstr(R,22),
                        'harmonic_strain_fraction':mp.nstr(harmonic_fraction,30),
                        'vortical_strain_fraction':mp.nstr(vortical_fraction,30),
                        'relative_error':mp.nstr(er,20)})

print(json.dumps({
 'precision_bits_requested':BITS,'mpmath_dps':mp.mp.dps,'status':'PASS',
 'interpretation':'The Hodge strain microscope is an exact physical-scale projector: the center harmonic strain is determined solely by the quadrupolar normal-flux moment on the sphere, while translation and higher harmonic boundary content are invisible. For an exact smooth divergence-free strain carrier tangent to its outer sphere, the microscope recovers S_h(R)=(1-R^2)S and S_v(R)=R^2 S, locating the transition from locally harmonic to locally vortical support without Fourier shells or an imposed source cutoff.',
 'projector_checks':projection_rows[:12],
 'carrier_profile':carrier
},indent=2))
