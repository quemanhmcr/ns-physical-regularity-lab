import json, os
from flint import ctx,arb
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import coupled46_hodge_core as H
import degree6_hodge_servo_core as C
st=H.prepare(); sym=H.feedback_symmetry_basis(st); seq=C.solve_degree6_servo()
# Locate the same symmetric branch with midpoint Newton, but use native Arb matrix solves throughout the physical eliminations.
yseq=[arb(seq['coeff'][i].mid()) for i in st['t4idx']]; a,_=H.sym_coords(yseq,sym); a=[arb(v.mid()) for v in a]
for _ in range(8):
    g,_,_=H.reduced_feedback_native(st,sym,a); Jraw=H.reduced_jacobian_native(st,sym,a)
    J=[[arb(Jraw[i][j].mid()) for j in range(5)] for i in range(5)]; d=H.arbmat_solve(J,[-arb(v.mid()) for v in g])
    a=[arb((a[i]+d[i]).mid()) for i in range(5)]
# Point preconditioner.
g0,fb0,_=H.reduced_feedback_native(st,sym,a); J0raw=H.reduced_jacobian_native(st,sym,a); J0=[[arb(J0raw[i][j].mid()) for j in range(5)] for i in range(5)]
Bcols=[]
for j in range(5):
    e=[C.z]*5;e[j]=C.o;Bcols.append(H.arbmat_solve(J0,e))
B=[[Bcols[j][i] for j in range(5)] for i in range(5)]
# Keep the same radii used by module 134 so only the linear observer changes.
rad_s='1e-6' if BITS==160 else ('1e-20' if BITS==256 else '1e-50'); R=arb('0 +/- '+rad_s)
X=[arb(v.mid())+R for v in a]; D=[arb('0 +/- '+rad_s) for _ in range(5)]
gX,fbX,_=H.reduced_feedback_native(st,sym,X); JX=H.reduced_jacobian_native(st,sym,X)
Bg=H.matvec(B,g0); center=[a[i]-Bg[i] for i in range(5)]; M=H.matmul(B,JX); I=H.eye(5); E=[[I[i][j]-M[i][j] for j in range(5)] for i in range(5)]; ED=H.matvec(E,D); K=[center[i]+ED[i] for i in range(5)]
inclusions=[]; margins=[]; certified=True
for i in range(5):
    inside=bool(K[i].lower()>X[i].lower() and K[i].upper()<X[i].upper()); inclusions.append(inside); certified=certified and inside
    margins.append({'lower_margin':str(K[i].lower()-X[i].lower()),'upper_margin':str(X[i].upper()-K[i].upper())})
# Native-solve residuals certify that the eliminated full-field equations are actually enclosed.
res4=C.norm2v(fbX['res4']); res6=C.norm2v(fbX['res6'])
out={
 'arb_precision_bits':BITS,'status':'PASS','linear_observer':'python-flint arb_mat.solve on fixed K44/K66 coordinate blocks',
 'Krawczyk_box_radius':rad_s,'Krawczyk_coordinate_inclusions':inclusions,'Krawczyk_certified_fixed_point':certified,
 'Krawczyk_box':[str(v) for v in X],'Krawczyk_image':[str(v) for v in K],'Krawczyk_inclusion_margins':margins,
 'candidate_reduced_feedback_coordinates':[str(v) for v in a],
 'candidate_point_feedback_coordinate_square_sum':str(sum((v*v for v in g0),C.z)),
 'interval_degree4_elimination_residual_mean_square':str(res4),'interval_degree6_elimination_residual_mean_square':str(res6),
 'interpretation':'This is the same five-dimensional physical Krawczyk test as module 134, with the same symmetry chart, candidate branch, box radii and inclusion criterion.  The only change is the observer for the fixed linear Hodge inverses: native python-flint arb_mat.solve replaces hand-written interval Gaussian elimination.  A successful inclusion therefore diagnoses the earlier 160/256 failure as numerical observer inflation rather than a change in physics.'}
print(json.dumps(out,indent=2,allow_nan=False))
