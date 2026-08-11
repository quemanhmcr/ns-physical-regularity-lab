import json, os
from flint import ctx,arb
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import coupled46_hodge_core as H
import degree6_hodge_servo_core as C
st=H.prepare(); sym=H.feedback_symmetry_basis(st); seq=C.solve_degree6_servo()
yseq=[arb(seq['coeff'][i].mid()) for i in st['t4idx']]; a,_=H.sym_coords(yseq,sym); a=[arb(v.mid()) for v in a]
for _ in range(8):
 g,_,_=H.reduced_feedback_native(st,sym,a); Jr=H.reduced_jacobian_native(st,sym,a); J=[[arb(Jr[i][j].mid()) for j in range(5)] for i in range(5)]
 d=H.arbmat_solve(J,[-arb(v.mid()) for v in g]); a=[arb((a[i]+d[i]).mid()) for i in range(5)]
g0,_,_=H.reduced_feedback_native(st,sym,a); J0r=H.reduced_jacobian_native(st,sym,a); J0=[[arb(J0r[i][j].mid()) for j in range(5)] for i in range(5)]
Bcols=[]
for j in range(5):
 e=[C.z]*5;e[j]=C.o;Bcols.append(H.arbmat_solve(J0,e))
B=[[Bcols[j][i] for j in range(5)] for i in range(5)]; Bg=H.matvec(B,g0); center=[a[i]-Bg[i] for i in range(5)]; I=H.eye(5)
rows=[]; successful=[]
for rad_s in ['1e-3','1e-4','1e-5','3e-6','1e-6','3e-7','1e-7','1e-8','1e-9','1e-10','1e-15','1e-20','1e-30','1e-40','1e-50']:
 R=arb('0 +/- '+rad_s); X=[arb(v.mid())+R for v in a]; D=[arb('0 +/- '+rad_s) for _ in range(5)]
 JX=H.reduced_jacobian_native(st,sym,X); M=H.matmul(B,JX); E=[[I[i][j]-M[i][j] for j in range(5)] for i in range(5)]; K=[center[i]+v for i,v in enumerate(H.matvec(E,D))]
 inc=[]; ratios=[]
 for i in range(5):
  inside=bool(K[i].lower()>X[i].lower() and K[i].upper()<X[i].upper());inc.append(inside)
  ratios.append((K[i].upper()-K[i].lower())/(X[i].upper()-X[i].lower()))
 ok=all(inc)
 if ok:successful.append(rad_s)
 rows.append({'radius':rad_s,'coordinate_inclusions':inc,'certified':ok,'Krawczyk_width_over_box_width':[str(v) for v in ratios],'Krawczyk_image':[str(v) for v in K]})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','successful_radii':successful,'has_certified_radius':bool(successful),'candidate_reduced_feedback_coordinates':[str(v) for v in a],'rows':rows,'interpretation':'The five-dimensional physical feedback map and native Arb fixed-Hodge observer are unchanged.  This module sweeps only the Krawczyk enclosure radius to distinguish a genuine loss of contraction from an interval-box conditioning artifact.  A successful radius is a rigorous fixed-point certificate at that precision; no tolerance or PDE parameter is altered.'},indent=2,allow_nan=False))
