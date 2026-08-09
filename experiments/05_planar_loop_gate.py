import json, os
from flint import arb, ctx

BITS = int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec = BITS
pi = arb.pi()

# Regularized Biot-Savart derivative evaluated by a periodic Riemann sum.
# For planar targets/sources the geometry forces u_theta to be normal to the plane,
# so X_theta dot u_theta must contain zero independently of quadrature resolution.

def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def dot(a,b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def add(a,b):
    return (a[0]+b[0],a[1]+b[1],a[2]+b[2])

def mul(c,a):
    return (c*a[0],c*a[1],c*a[2])

def sub(a,b):
    return (a[0]-b[0],a[1]-b[1],a[2]-b[2])

def ellipse(th,a,b):
    c,s=th.cos(),th.sin()
    return (a*c,b*s,arb(0)), (-a*s,b*c,arb(0))

def sigma_at(theta,a,b,Gamma,core,N):
    x, xt = ellipse(theta,a,b)
    h = 2*pi/N
    su = (arb(0),arb(0),arb(0))
    for j in range(N):
        ph = h*j
        xp, xpp = ellipse(ph,a,b)
        r = sub(x,xp)
        d2 = dot(r,r) + core*core
        d32 = d2*d2.sqrt()
        d52 = d32*d2
        c1 = cross(xpp,xt)
        c0 = cross(xpp,r)
        term = add(mul(1/d32,c1), mul(-3*dot(r,xt)/d52,c0))
        su = add(su, mul(h,term))
    uth = mul(Gamma/(4*pi),su)
    return dot(xt,uth)/dot(xt,xt)

aspects=['1','2','5','20','100']
core_fracs=['0.03','0.1','0.3']
scales=[('1','1'),('1e-6','1e9'),('1e6','1e-9')]
N=128
probes=12
rows=[]
for ar_s in aspects:
  ar=arb(ar_s)
  for cf_s in core_fracs:
    cf=arb(cf_s)
    for R_s,G_s in scales:
      R,G=arb(R_s),arb(G_s)
      a=R*ar.sqrt(); b=R/ar.sqrt()  # fixed area pi R^2 while varying anisotropy
      core=cf*R
      max_mid=0.0
      for k in range(probes):
        th = 2*pi*(arb(k)+arb('0.371'))/probes
        sig=sigma_at(th,a,b,G,core,N)
        if not sig.contains(0):
            raise AssertionError(f'planar stretching gate failed: aspect={ar_s}, core={cf_s}, R={R_s}, G={G_s}, k={k}, sigma={sig}')
        max_mid=max(max_mid,abs(float(sig.mid())))
      rows.append({'aspect':ar_s,'core_over_R':cf_s,'R':R_s,'Gamma':G_s,'max_sigma_mid_abs':max_mid})

print(json.dumps({
 'arb_precision_bits':BITS,
 'cases':len(rows),
 'probes_per_case':probes,
 'quadrature_nodes':N,
 'status':'PASS',
 'interpretation':'Every tested planar closed ellipse, including highly anisotropic ones, has zero tangential self-stretching under the regularized Biot-Savart geometry. Planar anisotropy alone is not the missing resource.',
 'rows':rows
},indent=2))
