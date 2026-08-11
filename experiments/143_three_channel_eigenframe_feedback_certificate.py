import json, os
from flint import arb,ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160:raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import coupled46_hodge_core as H
import degree6_hodge_servo_core as C
z=C.z;o=C.o;rt2=arb(2).sqrt()

def ppow(P,n):
 q={(0,0,0):o}
 for _ in range(n):q=C.pmul(q,P)
 return q
def substitute(P,L):
 out={};maxd=max((max(e) for e in P),default=0);cache=[[None]*(maxd+1) for _ in range(3)]
 for i in range(3):
  cache[i][0]={(0,0,0):o}
  for k in range(1,maxd+1):cache[i][k]=C.pmul(cache[i][k-1],L[i])
 for e,v in P.items():out=C.padd(out,C.pscale(v,C.pmul(C.pmul(cache[0][e[0]],cache[1][e[1]]),cache[2][e[2]])))
 return out
def group_action(V,R,X):
 L=[]
 for i in range(3):
  q={}
  for j in range(3):q=C.padd(q,C.pscale(R[i][j],X[j]))
  L.append(q)
 VR=[substitute(V[i],L) for i in range(3)];out=[]
 for i in range(3):
  q={}
  for j in range(3):q=C.padd(q,C.pscale(R[i][j],VR[j]))
  out.append(q)
 return tuple(out)
def mm(A,B):return tuple(tuple(sum((A[i][k]*B[k][j] for k in range(3)),z) for j in range(3)) for i in range(3))
def coords_in_basis(V,basis,d):
 cols=[C.flatten(v,d) for v in basis];sel,piv=C.independent(cols)
 if len(sel)!=len(basis):raise AssertionError(('basis rank',len(sel),len(basis)))
 A=[[cols[j][piv[i]] for j in range(len(basis))] for i in range(len(basis))];b=C.flatten(V,d);c=H.arbmat_solve(A,[b[p] for p in piv]);res=C.vadd(V,C.vscale(-1,C.combine(c,basis)))
 return c,C.norm2v(res)

st=H.prepare();Y=[st['V6b'][i] for i in st['t4idx']]
R1=((z,-o,z),(-o,z,z),(z,z,-o));t=rt2-o;den=o+t*t;v=(o/rt2,o/rt2,t);R2=tuple(tuple(arb(2)*v[i]*v[j]/den-(o if i==j else z) for j in range(3)) for i in range(3));R3=mm(R1,R2)
# Group-average the complete 9D T4 feedback basis; its image is the full D2-fixed T4 subspace.
avg=[]
for V in Y:
 q=V
 for R in (R1,R2,R3):q=C.vadd(q,group_action(V,R,st['X']))
 avg.append(C.vscale(arb(1)/4,q))
sel,_=C.independent([C.flatten(v,6) for v in avg]);B3=[avg[i] for i in sel]
if len(B3)!=3:raise AssertionError(('D2 fixed T4 dimension',len(B3)))
# Express the three physical fields in the original nine T4 coordinates.
B3coords=[]
for V in B3:
 c,res=coords_in_basis(V,Y,6)
 if not res.contains(0):raise AssertionError(('B3 coordinate residual',res))
 B3coords.append(c)
# Coordinate extraction in B3.
cols3=[C.flatten(v,6) for v in B3];_,piv3=C.independent(cols3);A3=[[cols3[j][piv3[i]] for j in range(3)] for i in range(3)]
def y_from_a(a):
 return [sum((a[k]*B3coords[k][i] for k in range(3)),z) for i in range(9)]
def a_from_y(y):
 V=C.combine(y,Y);b=C.flatten(V,6);a=H.arbmat_solve(A3,[b[p] for p in piv3]);res=C.vadd(V,C.vscale(-1,C.combine(a,B3)));return a,res
def G(a):
 y=y_from_a(a);fb=H.feedback_map_native(st,y);g,res=a_from_y(fb['F']);return g,fb,res
def JG(a):
 y=y_from_a(a);J9=H.feedback_jacobian_native(st,y);cols=[]
 for k in range(3):
  dy=B3coords[k];df=[sum((J9[i][j]*dy[j] for j in range(9)),z) for i in range(9)];c,res=a_from_y(df)
  if not C.norm2v(res).contains(0):raise AssertionError(('D2 Jacobian closure',k,C.norm2v(res)))
  cols.append(c)
 return [[cols[j][i] for j in range(3)] for i in range(3)]
# Start from the already certified 5D branch midpoint and project to D2 coordinates.
sym=H.feedback_symmetry_basis(st);seq=C.solve_degree6_servo();ys=[arb(seq['coeff'][i].mid()) for i in st['t4idx']];a5,_=H.sym_coords(ys,sym);a5=[arb(v.mid()) for v in a5]
for _ in range(8):
 g,_,_=H.reduced_feedback_native(st,sym,a5);Jr=H.reduced_jacobian_native(st,sym,a5);J=[[arb(Jr[i][j].mid()) for j in range(5)] for i in range(5)];dd=H.arbmat_solve(J,[-arb(v.mid()) for v in g]);a5=[arb((a5[i]+dd[i]).mid()) for i in range(5)]
y5=H.y_from_sym(a5,sym);a,res=a_from_y(y5)
if not C.norm2v(res).contains(0):raise AssertionError(('certified branch not D2 fixed',C.norm2v(res)))
a=[arb(v.mid()) for v in a]
for _ in range(6):
 g,_,_=G(a);Jr=JG(a);J=[[arb(Jr[i][j].mid()) for j in range(3)] for i in range(3)];dd=H.arbmat_solve(J,[-arb(v.mid()) for v in g]);a=[arb((a[i]+dd[i]).mid()) for i in range(3)]
# Krawczyk in the three physical eigenframe feedback channels.
g0,fb0,_=G(a);J0r=JG(a);J0=[[arb(J0r[i][j].mid()) for j in range(3)] for i in range(3)];Bcols=[]
for j in range(3):e=[z]*3;e[j]=o;Bcols.append(H.arbmat_solve(J0,e))
Binv=[[Bcols[j][i] for j in range(3)] for i in range(3)];rad_s='1e-20';R=arb('0 +/- '+rad_s);X=[arb(v.mid())+R for v in a];D=[arb('0 +/- '+rad_s) for _ in range(3)];JX=JG(X);Bg=H.matvec(Binv,g0);center=[a[i]-Bg[i] for i in range(3)];M=H.matmul(Binv,JX);I=H.eye(3);E=[[I[i][j]-M[i][j] for j in range(3)] for i in range(3)];K=[center[i]+v for i,v in enumerate(H.matvec(E,D))];inc=[bool(K[i].lower()>X[i].lower() and K[i].upper()<X[i].upper()) for i in range(3)]
if not all(inc):raise AssertionError(('3D D2 Krawczyk',inc,K,X))
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','full_T4_feedback_dimension':9,'D2_fixed_T4_feedback_dimension':3,'candidate_distance_from_D2_fixed_subspace_mean_square':str(C.norm2v(res)),'D2_feedback_coordinates':[str(v) for v in a],'Krawczyk_box_radius':rad_s,'Krawczyk_coordinate_inclusions':inc,'Krawczyk_certified_fixed_point':True,'Krawczyk_image':[str(v) for v in K],'interpretation':'The complete nine-dimensional degree-six T4 feedback sector is group-averaged under the proper D2 eigenframe stabilizer of the stationary strain.  Its physical fixed subspace has dimension three, exactly the even-l formula (l+2)/2 at l=4.  The previously certified coupled branch lies in this three-dimensional space.  Re-expressing the exact Hodge/Euler feedback map in these three physical eigenframe channels and applying Krawczyk again certifies the same local coupled servo.  Thus the apparent feedback complexity reduces 9 -> 5 -> 3 when the full physical stabilizer is respected.'},indent=2,allow_nan=False))
