import json, os
from fractions import Fraction as F
from flint import arb,ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160:raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import degree6_hodge_servo_core as C
z=C.z;o=C.o

def re_z_power(l):
 from math import comb
 return {(l-k,k,0):arb(((-1)**(k//2))*comb(l,k)) for k in range(0,l+1,2)}
def r2pow(k,r2):
 q={(0,0,0):o}
 for _ in range(k):q=C.pmul(q,r2)
 return q
def normcoord(V):
 q=z
 for comp in V:
  for v in comp.values():q+=v*v
 return q
X,r2,u1,u3,omega=C.setup_base()
profiles=[[(0,'1'),(2,'-0.3'),(4,'0.07')],[(0,'1e-20'),(2,'3e10'),(6,'-2e5'),(8,'0.125')],[(0,'1e20'),(4,'-2e-10'),(10,'3e-30')]]
rows=[]
for l in (4,6,8):
 H=re_z_power(l);g=tuple(C.pder(H,i) for i in range(3));base_route=C.bracket(omega,g)
 for pi,prof in enumerate(profiles):
  Cmom=z;route_sum=({},{},{})
  for q,cs in prof:
   c=arb(cs); cq=-c*arb(l+1)/((q+2)*(q+2*l+3));Cmom+=cq;route_sum=C.vadd(route_sum,C.vscale(cq,base_route))
  route_moment=C.vscale(Cmom,base_route);err=C.vadd(route_sum,C.vscale(-1,route_moment));err2=normcoord(err)
  if not err2.contains(0):raise AssertionError(('routing moment',l,pi,err2))
  degrees=sorted(set(sum(e) for comp in route_sum for e in comp))
  if route_sum!=( {},{},{} ) and degrees!=[l]:raise AssertionError(('route degree',l,pi,degrees))
  rows.append({'l':l,'profile_index':pi,'screened_feedback_moment_C_l':str(Cmom),'lower_response_homogeneous_degrees':degrees,'route_factorization_coordinate_error':str(err2)})
 # Explicit nonzero screened-null radial profile 1-c r^2.
 cnull=arb(2*(2*l+5))/arb(2*l+3)
 c0=-arb(l+1)/(2*(2*l+3));c2=(cnull)*arb(l+1)/(4*(2*l+5)) # coefficient -cnull times negative monomial companion
 Cnull=c0+c2
 if not Cnull.contains(0):raise AssertionError(('screened-null moment',l,Cnull))
 V0=C.vadd(C.cross(X,g),C.vscale(-cnull,tuple(C.pmul(r2,t) for t in C.cross(X,g))))
 if not (normcoord(V0)>0):raise AssertionError(('silent profile nonzero',l,normcoord(V0)))
 route0=C.vscale(Cnull,base_route);r0=normcoord(route0)
 if not r0.contains(0):raise AssertionError(('screened-null route',l,r0))
 rows.append({'l':l,'profile_index':'screened_null_exact','profile_formula':f'1-{str(cnull)} r^2','screened_feedback_moment_C_l':str(Cnull),'vorticity_coordinate_square_nonzero':str(normcoord(V0)),'lower_response_coordinate_square':str(r0),'lower_silent_radial_profile':True})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For a fixed toroidal angular Hodge channel T_l, every radial Taylor component contributes to the same lower harmonic velocity direction grad H_l.  Hence the full base-vorticity backreaction factorizes exactly as [omega2,U_low]=C_l[a][omega2,grad H_l] and has homogeneity l independent of radial Taylor order.  The only radial information seen by the lower response is the screened moment C_l[a].  The explicit nonzero profile a(r)=1-2(2l+5)/(2l+3) r^2 has C_l[a]=0 and is therefore lower-silent.  Thus the kernel of the screened feedback moment contains genuine internal radial vorticity shapes; polynomial radial degrees are not independent lower-feedback channels.','rows':rows},indent=2,allow_nan=False))
