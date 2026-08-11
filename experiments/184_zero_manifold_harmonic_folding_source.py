import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import degree6_hodge_servo_core as C
z=C.z;one=C.o
X=({(1,0,0):one},{(0,1,0):one},{(0,0,1):one});x,y,zz=X

def eval_at_axis(P,s):
 out=z
 for (a,b,c),v in P.items():
  if a==0 and b==0:out+=v*(s**c)
 return out
def vaxis(V,s):return tuple(eval_at_axis(q,s) for q in V)
def d2z(V):return tuple(C.pder(C.pder(q,2),2) for q in V)
def vn2(v):return sum((a*a for a in v),z)

def harmonic_phi_m(m):
 # phi=-Re[(x+i z)^(m+1)]/(m+1), m even. On z-axis u_x=(-1)^(m/2+1) z^m, u_z=0.
 from math import comb
 n=m+1;P={}
 for k in range(0,n+1,2):
  # term x^(n-k) (i z)^k
  P[(n-k,0,k)]=arb(-(((-1)**(k//2))*comb(n,k)))/arb(n)
 return P
rows=[]
# Affine harmonic potential phi=xz: Hessian of velocity is zero, no curvature creation from a straight line.
phi_aff=C.pmul(x,zz);u_aff=tuple(C.pder(phi_aff,i) for i in range(3));curl_aff=C.curl(u_aff);div_aff=C.div(u_aff)
if not C.savg(C.vdot(curl_aff,curl_aff)).contains(0) or not C.savg(C.pmul(div_aff,div_aff)).contains(0):raise AssertionError('affine potential')
Kaff=d2z(u_aff)
for q in Kaff:
 if q: raise AssertionError(('affine K source',Kaff))
for m in (2,4,6,8,12):
 phi=harmonic_phi_m(m);lap=C.plap(phi)
 if any(not v.contains(0) for v in lap.values()):raise AssertionError(('phi harmonic',m,lap))
 u=tuple(C.pder(phi,i) for i in range(3));curl=C.curl(u);div=C.div(u)
 if not C.savg(C.vdot(curl,curl)).contains(0):raise AssertionError(('curl',m))
 if not C.savg(C.pmul(div,div)).contains(0):raise AssertionError(('div',m))
 Ksrc=d2z(u)
 for ss in ('0.25','0.5','1'):
  s=arb(ss);ua=vaxis(u,s);ka=vaxis(Ksrc,s)
  expected=arb(m*(m-1))*(s**(m-2))
  if not (abs(ka[0])/expected).contains(1):raise AssertionError(('curvature source',m,ss,ka,expected))
  if not ka[1].contains(0) or not ka[2].contains(0):raise AssertionError(('transverse source',m,ss,ka))
  rows.append({'harmonic_velocity_axis_power_m':m,'axis_parameter_z':ss,'axis_velocity':[str(a) for a in ua],'straight_line_curvature_source_DtK_at_t0':[str(a) for a in ka],'curvature_source_abs_expected_m_mminus1_z_pow_mminus2':str(expected),'vorticity_mean_square':'0','velocity_divergence_mean_square':'0'})
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'affine_harmonic_straight_line_curvature_source':'0',
 'interpretation':(
  'For any material curve X(s,t), exact differentiation gives D_t T=(grad u)T and D_t K=(grad u)K+(grad grad u)[T,T].  Thus a straight zero line K=0 can be folded only by the velocity Hessian; a common affine strain has zero curvature source. '
  'The artifact then constructs explicit non-affine harmonic potentials phi=-Re[(x+i z)^(m+1)]/(m+1).  Their velocity is incompressible and irrotational, yet on the initially straight z-axis the transverse velocity is proportional to z^m and the instantaneous curvature source is exactly m(m-1)z^(m-2).  Hence degenerate vorticity-zero manifolds can be folded by an irrotational Hodge-harmonic actor without any local vorticity or Kelvin mutation. '
  'This kills any claim that fixed-time multiplicity must itself pay a null-vorticity ancestry cost.  The folding actor is instead a non-affine harmonic velocity field.  Because a nontrivial global harmonic finite-energy field cannot persist on all space, a physical realization must be sourced from outside the local ball; the next attack must measure its Hodge occupancy/source radius rather than assign it an enstrophy tax.'),
 'rows':rows
},indent=2,allow_nan=False))
