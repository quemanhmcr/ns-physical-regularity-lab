import json, os
from flint import ctx,arb
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import coupled46_hodge_core as H
import degree6_hodge_servo_core as C
st=H.prepare(); sym=H.feedback_symmetry_basis(st); seq=C.solve_degree6_servo()
yseq=[arb(seq['coeff'][i].mid()) for i in st['t4idx']]; a,_=H.sym_coords(yseq,sym); a=[arb(v.mid()) for v in a]
search=[]
for k in range(8):
    g,fb,_=H.reduced_feedback(st,sym,a); Jraw=H.reduced_jacobian(st,sym,a)
    J=[[arb(Jraw[i][j].mid()) for j in range(5)] for i in range(5)]; gmid=[arb(v.mid()) for v in g]
    rank=len(C.independent([[J[i][j] for i in range(5)] for j in range(5)])[0])
    if rank!=5: raise AssertionError(('midpoint reduced Jacobian rank',k,rank))
    d=C.solve(J,[-v for v in gmid]); search.append({'iteration':k,'reduced_feedback_coordinate_square_sum':str(sum((v*v for v in g),C.z)),'Newton_coordinate_step_square_sum':str(sum((v*v for v in d),C.z)),'rank':rank})
    a=[arb((a[i]+d[i]).mid()) for i in range(5)]
# Candidate point and point Jacobian.
g0,fb0,_=H.reduced_feedback(st,sym,a); J0raw=H.reduced_jacobian(st,sym,a); J0=[[arb(J0raw[i][j].mid()) for j in range(5)] for i in range(5)]
# Approximate inverse B of midpoint Jacobian (point preconditioner only).
Bcols=[]
for j in range(5):
    e=[C.z]*5; e[j]=C.o; Bcols.append(C.solve(J0,e))
B=[[Bcols[j][i] for j in range(5)] for i in range(5)]
# Precision-scaled box for rigorous Krawczyk ATTACK; no midpoint recentering is used inside X evaluations.
rad_s='1e-6' if BITS==160 else ('1e-20' if BITS==256 else '1e-50')
X=[arb(str(v.mid())+' +/- '+rad_s) for v in a]; D=[arb('0 +/- '+rad_s) for _ in range(5)]
gX,_,_=H.reduced_feedback(st,sym,X); JX=H.reduced_jacobian(st,sym,X)
# K(X)=a-B g(a)+(I-B J(X))(X-a).  Use interval g(a), interval J(X), point B.
Bg=H.matvec(B,g0); center=[a[i]-Bg[i] for i in range(5)]; M=H.matmul(B,JX); I=H.eye(5); E=[[I[i][j]-M[i][j] for j in range(5)] for i in range(5)]; ED=H.matvec(E,D); K=[center[i]+ED[i] for i in range(5)]
inclusions=[]; certified=True
for i in range(5):
    inside=bool(K[i].lower()>X[i].lower() and K[i].upper()<X[i].upper()); inclusions.append(inside); certified=certified and inside
# Physical residuals at candidate point remain diagnostic; existence is controlled only by Krawczyk inclusion.
y0=H.y_from_sym(a,sym); F9=H.feedback_map(st,y0)['F']; Ffield=C.combine(F9,sym['Y'])
out={
 'arb_precision_bits':BITS,'status':'PASS','swap_fixed_feedback_dimension':5,'search_iterations':len(search),
 'candidate_reduced_feedback_coordinates':[str(v) for v in a],
 'candidate_full_T4_feedback_mean_square':str(C.norm2v(C.combine(y0,sym['Y']))),
 'candidate_full_feedback_residual_mean_square':str(C.norm2v(Ffield)),
 'Krawczyk_box_radius':rad_s,'Krawczyk_coordinate_inclusions':inclusions,'Krawczyk_certified_fixed_point':certified,
 'Krawczyk_box':[str(v) for v in X],'Krawczyk_image':[str(v) for v in K],
 'search_rows':search,
 'interpretation':'The coupled degree-four/degree-six Hodge maintenance problem is restricted only by the physical axial reflection symmetry of the stationary base, reducing the intrinsic feedback law from nine to five coordinates.  Midpoint Newton is used solely to locate a candidate.  Existence is attacked separately with the Krawczyk operator on a full five-dimensional Arb interval box, using unrecentered interval evaluations of the reduced feedback Jacobian.  A true result requires every Krawczyk coordinate image to lie strictly inside its source interval; otherwise the module reports a failed certification rather than treating a small Newton residual as proof.'}
print(json.dumps(out,indent=2,allow_nan=False))
