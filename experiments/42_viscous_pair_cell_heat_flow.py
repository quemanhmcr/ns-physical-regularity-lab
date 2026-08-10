import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def det(a,b,c): return dot(a,cross(b,c))
def vadd(a,b): return tuple(a[i]+b[i] for i in range(3))
def vscale(c,a): return tuple(c*a[i] for i in range(3))
def certify_one(x,label):
    tol=arb('1e-30')
    if not x.contains(1) or not (x > 1-tol and x < 1+tol):
        raise AssertionError((label,'ratio not tightly certified around one',x))

# Exact NS solution u(y,t)=(A e^{-nu k^2 t} cos ky,0,B e^{-nu m^2 t} cos my).
# Nonlinearity vanishes because u_y=0 and all fields depend only on y.
k=arb(1); m=arb(2)
y1=pi/6; y2=pi/3; dy=y2-y1
As=['1e-30','1','1e30']
Bs=['1e-30','1','1e30']
nus=['1e-30','1','1e30']
qs=['0','1e-12','0.01','0.1','1','10','100']  # q=nu t
rows=[]
for As_ in As:
  A=arb(As_)
  for Bs_ in Bs:
    B=arb(Bs_)
    # D0 at q=0, using R_y only; R_x,R_z never enter this determinant.
    def omega(y,q):
      Ek=(-(k*k)*q).exp(); Em=(-(m*m)*q).exp()
      return (-B*m*Em*(m*y).sin(),arb(0),A*k*Ek*(k*y).sin())
    p0=omega(y1,arb(0)); r0=omega(y2,arb(0))
    R=(arb('0.37'),dy,arb('-0.21'))
    D0=det(p0,R,r0)
    if D0.contains(0): raise AssertionError(('heat-flow base D degenerate',As_,Bs_,D0))
    for nus_ in nus:
      nu=arb(nus_)
      for qs_ in qs:
        q=arb(qs_)
        p=omega(y1,q); rr=omega(y2,q)
        D=det(p,R,rr)
        expected_ratio=(-(k*k+m*m)*q).exp()
        ratio=(D/D0)/expected_ratio
        certify_one(ratio,('viscous pair-cell decay',As_,Bs_,nus_,qs_))
        # Vortex stretching and bridge inhomogeneity are exactly zero: grad u has only column y, omega_y=0.
        # Laplacians act with eigenvalues -m^2 on omega_x and -k^2 on omega_z.
        lap_p=(-m*m*p[0],arb(0),-k*k*p[2])
        lap_r=(-m*m*rr[0],arb(0),-k*k*rr[2])
        visc=nu*(det(lap_p,R,rr)+det(p,R,lap_r))
        expected_rate=-nu*(k*k+m*m)
        rate=visc/D
        certify_one(rate/expected_rate,('viscous Ddot/D rate',As_,Bs_,nus_,qs_))
        # Direct material derivative uses the actual differential advection of the exact heat flow.
        Ek=(-(k*k)*q).exp(); Em=(-(m*m)*q).exp()
        Rdot=(A*Ek*((k*y2).cos()-(k*y1).cos()),arb(0),B*Em*((m*y2).cos()-(m*y1).cos()))
        pdot=vscale(nu,lap_p); rdot=vscale(nu,lap_r)
        direct=det(pdot,R,rr)+det(p,Rdot,rr)+det(p,R,rdot)
        certify_one(direct/visc,('direct heat pair balance',As_,Bs_,nus_,qs_))
        rows.append({
          'A':As_,'B':Bs_,'nu':nus_,'q_nu_t':qs_,
          'D_over_D0':str(D/D0),'expected_exp_minus_5q':str(expected_ratio),
          'Ddot_over_D':str(rate),'expected_minus_5nu':str(expected_rate),
          'bridge_current':'0','direct_over_viscous_current':str(direct/visc),
        })

print(json.dumps({
 'arb_precision_bits':BITS,
 'status':'PASS',
 'cases':len(rows),
 'interpretation':(
   'An exact nonlinear-null Navier-Stokes heat flow isolates the viscous branch of the material pair-cell law. '
   'Because grad u has only a y column while omega_y=0, vortex stretching and bridge-gradient renewal vanish identically. '
   'The noncoplanar pair ancestry cell decays exactly as exp[-nu(k^2+m^2)t]=exp[-5 nu t], and Ddot/D=-5 nu, certifying the viscous pair-cell current independently of the nonlinear bridge branch.'
 ),
 'rows':rows,
},indent=2,allow_nan=False))
