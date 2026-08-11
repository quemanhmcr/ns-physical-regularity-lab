import json, os
from flint import ctx,arb
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import coupled46_hodge_core as H
import degree6_hodge_servo_core as C
st=H.prepare();sym=H.feedback_symmetry_basis(st);seq=C.solve_degree6_servo(); y=[arb(seq['coeff'][i].mid()) for i in st['t4idx']];a,_=H.sym_coords(y,sym);a=[arb(v.mid()) for v in a]
for _ in range(8):
 g,_,_=H.reduced_feedback_native(st,sym,a);Jr=H.reduced_jacobian_native(st,sym,a);J=[[arb(Jr[i][j].mid()) for j in range(5)] for i in range(5)];d=H.arbmat_solve(J,[-arb(v.mid()) for v in g]);a=[arb((a[i]+d[i]).mid()) for i in range(5)]
# Common certified radius from module136 at all 160/256/512 precisions.
rad_s='1e-20';R=arb('0 +/- '+rad_s);X=[arb(v.mid())+R for v in a];D=[arb('0 +/- '+rad_s) for _ in range(5)]
g0,_,_=H.reduced_feedback_native(st,sym,a);J0r=H.reduced_jacobian_native(st,sym,a);J0=[[arb(J0r[i][j].mid()) for j in range(5)] for i in range(5)]
Bcols=[]
for j in range(5):e=[C.z]*5;e[j]=C.o;Bcols.append(H.arbmat_solve(J0,e))
B=[[Bcols[j][i] for j in range(5)] for i in range(5)];JX=H.reduced_jacobian_native(st,sym,X);Bg=H.matvec(B,g0);center=[a[i]-Bg[i] for i in range(5)];M=H.matmul(B,JX);I=H.eye(5);E=[[I[i][j]-M[i][j] for j in range(5)] for i in range(5)];K=[center[i]+v for i,v in enumerate(H.matvec(E,D))]
inc=[bool(K[i].lower()>X[i].lower() and K[i].upper()<X[i].upper()) for i in range(5)]
if not all(inc):raise AssertionError(('common root box not Krawczyk certified',BITS,inc,K,X))
# Evaluate full higher response on the entire root-containing box.  Positive lower norm then holds at the actual unique root.
_,fbX,_=H.reduced_feedback_native(st,sym,X); Hresp=H.higher_responses_from_coupled(st,fbX);rows=[]
for d in (8,10,12):
 RR,PP,NN=Hresp[d];tot=C.norm2v(RR);prod=C.norm2v(PP);null=C.norm2v(NN)
 if not (null>0):raise AssertionError(('higher null response not separated from zero on certified root box',d,null))
 rows.append({'degree':d,'total_response_mean_square':str(tot),'productive_projection_mean_square':str(prod),'transaction_null_mean_square':str(null),'transaction_null_strictly_positive_on_entire_certified_root_box':True})
out={'arb_precision_bits':BITS,'status':'PASS','certified_root_box_radius':rad_s,'Krawczyk_inclusions':inc,'candidate_coordinates':[str(v) for v in a],'degree4_interval_elimination_residual_mean_square':str(C.norm2v(fbX['res4'])),'degree6_interval_elimination_residual_mean_square':str(C.norm2v(fbX['res6'])),'interpretation':'The unique swap-symmetric coupled degree-four/degree-six Hodge servo is enclosed in a common radius 1e-20 Krawczyk box that is independently recertified at each precision.  The complete Euler vorticity interaction of omega2+V4+V6 with its exact Hodge velocity is then grouped by physical homogeneity.  Degrees four and six are the maintained levels.  The remaining degree-eight, ten and twelve responses are projected through the same sharp transaction projector.  A strictly positive lower bound for each null norm on the entire certified root box means every actual coupled root in that box necessarily emits those higher transaction-null responses; this is stronger than evaluating a midpoint candidate.','rows':rows}
print(json.dumps(out,indent=2,allow_nan=False))
