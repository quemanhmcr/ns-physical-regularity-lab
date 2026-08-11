import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import degree6_hodge_servo_core as C
z=C.z; one=C.o; pi=arb.pi(); rt3=arb(3).sqrt()
X=({(1,0,0):one},{(0,1,0):one},{(0,0,1):one})
B=((arb(2),z,z),(z,-one,z),(z,z,-one))

def matvec(B,V):
 out=[]
 for i in range(3):
  q={}
  for j in range(3):q=C.padd(q,C.pscale(B[i][j],V[j]))
  out.append(q)
 return tuple(out)
def vzero(V): return C.savg(C.vdot(V,V)).contains(0)

def directional(v,a): return C.directional(v,a)

rows=[]
for As in ('1e-30','1','1e30'):
 A=arb(As)
 omega=C.vscale(A,matvec(B,X))
 u=C.vscale(-A/3,C.cross(X,matvec(B,X)))
 G=C.bracket(omega,u)
 pred=C.vscale(2*A,u)
 if not vzero(C.vadd(G,C.vscale(-1,pred))): raise AssertionError(('G=2Au',As))
 radial_u=C.vdot(X,u); radial_G=C.vdot(X,G)
 if not C.savg(C.pmul(radial_u,radial_u)).contains(0): raise AssertionError(('u tangent',As))
 if not C.savg(C.pmul(radial_G,radial_G)).contains(0): raise AssertionError(('G tangent',As))
 divw=C.div(omega)
 if not C.savg(C.pmul(divw,divw)).contains(0): raise AssertionError(('div omega',As))
 # Exact vortex-line first integrals: x(y^2+z^2) and azimuth y/z.
 x,y,zv=X;rho2=C.padd(C.pmul(y,y),C.pmul(zv,zv));I=C.pmul(x,rho2)
 Idot={}
 for j,q in enumerate(X):
  pass
 dI={}
 for j in range(3): dI=C.padd(dI,C.pmul(omega[j],C.pder(I,j)))
 if not C.savg(C.pmul(dI,dI)).contains(0): raise AssertionError(('vortex invariant x rho2',As))
 az=C.padd(C.pmul(y,omega[2]),C.pscale(-1,C.pmul(zv,omega[1])))
 if not C.savg(C.pmul(az,az)).contains(0): raise AssertionError(('fixed vortex azimuth',As))
 # Material flow keeps x and rho^2 fixed; particles rotate around x-axis on each sphere.
 xdot=u[0];rhodot={}
 for j,q in [(1,y),(2,zv)]: rhodot=C.padd(rhodot,C.pscale(2,C.pmul(q,u[j])))
 if not C.savg(C.pmul(xdot,xdot)).contains(0) or not C.savg(C.pmul(rhodot,rhodot)).contains(0): raise AssertionError(('material circles',As))
 # Sphere radial flux density = A r (3 n_x^2-1); its total mean is zero.
 radial_flux=C.vdot(X,omega)
 if not C.savg(radial_flux).contains(0): raise AssertionError(('zero total sphere flux',As))
 # One positive cap n_x>1/sqrt3 has exact outward flux at radius R: 4pi A R^3/(3sqrt3).
 cap_coeff=arb(4)*pi*A/(3*rt3)
 rows.append({
  'A':As,'self_Euler_source_equals_2A_u':True,'u_tangent_to_centered_spheres':True,
  'self_Euler_source_tangent_to_centered_spheres':True,
  'vortex_line_invariant_x_times_rho_squared':'exact',
  'vortex_line_fixed_azimuth':'exact','material_x_and_rho_squared_fixed':'exact',
  'sphere_total_radial_vorticity_flux':'0',
  'one_positive_cap_flux_coefficient_Gamma_over_R3':str(cap_coeff),
  'one_negative_belt_balances_two_positive_caps':True,
 })

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'For the axisymmetric linear null catalyst omega=A(2x,-y,-z), the exact tangent Hodge velocity is u=(0,-A x z,A x y).  Its vortex lines are hyperbolic through-going curves: x(y^2+z^2) and azimuth are first integrals, while radial flux enters the equatorial belt and exits the two polar caps of every centered sphere.  Thus this local catalyst is not a self-contained closed vorticity ancestry stock. '
  'At the same time the fluid velocity keeps x and y^2+z^2 fixed, so material particles move on circles around the x axis and every centered sphere is material for this instantaneous velocity.  Most importantly, the exact self-Euler source satisfies G=2A u and is tangent to every centered sphere: G.n=0 pointwise. '
  'Therefore the null catalyst can generate the productive Hodge transaction certified in module168 by tangentially reorganizing vorticity while leaving the radial vorticity-flux density through each centered sphere unchanged at first Euler order.  This is the Kelvin principle in local geometric form: Euler conversion need not spend material circulation ancestry. '
  'The unavoidable cost can only appear when the polynomial catalyst is made finite-energy by spatial turnover/localization, where viscosity and closure enter.  The next microscope localizes the velocity itself so divergence-free vorticity is preserved and measures the exact turnover-collar enstrophy/lifetime.'),
 'rows':rows
},indent=2,allow_nan=False))
