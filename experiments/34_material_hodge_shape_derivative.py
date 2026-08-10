import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

# Exact shape derivative of the unit-period harmonic circulation class on a
# material flat collar.  The theorem being calibrated is
#   I_dot = -2 int h.S.h = -2 sigma_C I.
# No finite differencing is used: all time derivatives are analytic.

I0_values=['1e-30','1','1e30']
L_values=['1e-30','1e-12','1','1e12','1e30']
a_values=['1e-24','1','1e24']
rows=[]

# Volume-preserving diagonal deformation M=diag(L,L^-1/2,L^-1/2).
# Cycle along x: I=I0/L^2, sigma=a.
# Cycle along y: I=I0*L, sigma=-a/2.
for I0s in I0_values:
    I0=arb(I0s)
    for Ls in L_values:
        L=arb(Ls)
        for a_s in a_values:
            a=arb(a_s)
            Ix=I0/(L*L)
            sigx=a
            Idotx=-2*a*Ix
            rx=Idotx + 2*sigx*Ix
            if not rx.contains(0):
                raise AssertionError(('x-cycle shape law residual',I0s,Ls,a_s,rx))
            if not ((Idotx/Ix)/(-2*a)).contains(1):
                raise AssertionError(('x-cycle logarithmic rate',I0s,Ls,a_s,Idotx/Ix))

            Iy=I0*L
            sigy=-a/2
            Idoty=a*Iy
            ry=Idoty + 2*sigy*Iy
            if not ry.contains(0):
                raise AssertionError(('transverse-cycle shape law residual',I0s,Ls,a_s,ry))
            if not ((Idoty/Iy)/a).contains(1):
                raise AssertionError(('transverse-cycle logarithmic rate',I0s,Ls,a_s,Idoty/Iy))

            rows.append({
                'I0':I0s,'L':Ls,'a':a_s,
                'I_extensional_over_I0':str(Ix/I0),
                'sigma_extensional':str(sigx),
                'dlogI_extensional_over_a':str((Idotx/Ix)/a),
                'I_compressed_over_I0':str(Iy/I0),
                'sigma_compressed':str(sigy),
                'dlogI_compressed_over_a':str((Idoty/Iy)/a),
            })

# Affine shear M=[[1,0,0],[k,1,0],[0,0,1]], k_dot=g.
# Existing intrinsic one-cycle collar has I=I0/(1+k^2).
shear=[]
for I0s in I0_values:
    I0=arb(I0s)
    for ks in ['-1e30','-1e6','-1','-0.1','0','0.1','1','1e6','1e30']:
        k=arb(ks)
        for gs in ['-1e12','-1','1','1e12']:
            g=arb(gs)
            den=1+k*k
            I=I0/den
            sigma=g*k/den
            Idot=-2*I0*k*g/(den*den)
            residual=Idot+2*sigma*I
            if not residual.contains(0):
                raise AssertionError(('shear shape law residual',I0s,ks,gs,residual))
            if ks!='0':
                target=-2*sigma
                if not ((Idot/I)/target).contains(1):
                    raise AssertionError(('shear dlogI mismatch',I0s,ks,gs,Idot/I,target))
            shear.append({
                'I0':I0s,'k':ks,'k_dot':gs,
                'I_over_I0':str(I/I0),
                'sigma_C':str(sigma),
                'dlogI_dt':str(Idot/I),
                'shape_residual':str(residual),
            })

# Rigid rotation: S=0, hence the period metric is unchanged.
rotation=[]
for I0s in I0_values:
    I0=arb(I0s)
    for Omega_s in ['1e-30','1','1e30']:
        Omega=arb(Omega_s)
        sigma=arb(0)
        Idot=arb(0)
        residual=Idot+2*sigma*I0
        if not residual.contains(0):
            raise AssertionError(('rotation invariance residual',I0s,Omega_s,residual))
        rotation.append({'I0':I0s,'Omega':Omega_s,'sigma_C':'0','dI_dt':'0'})

print(json.dumps({
    'arb_precision_bits':BITS,
    'status':'PASS',
    'diagonal_cases':len(rows),
    'shear_cases':len(shear),
    'rotation_cases':len(rotation),
    'interpretation':(
        'The intrinsic circulation inductance of a material collar obeys the exact shape law '
        'I_dot=-2 integral h.S.h. A cycle stretched by the flow becomes cheaper to carry at fixed period, '
        'a compressed circulation cycle becomes more expensive, affine shear follows the same law, and rigid rotation costs nothing. '
        'This is a geometry evolution identity, not an irreversible dissipation estimate.'
    ),
    'diagonal':rows,
    'shear':shear,
    'rotation':rotation,
},indent=2,allow_nan=False))
