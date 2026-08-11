import json, os
from flint import ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import coupled46_hodge_core as H
import degree6_hodge_servo_core as C
st=H.prepare()
# The entire lower backreaction is carried by the degree-six T4 sector.
silent_back=[st['K46'][i] for i in st['silent']]
max_silent=C.z
for v in silent_back:
    q=C.norm2v(v)
    if not q.contains(0): raise AssertionError(('lower-silent control backreacts',q))
    if q>max_silent:max_silent=q
cols46=[C.flatten(st['K46'][i],4) for i in st['t4idx']]
r46=len(C.independent(cols46)[0])
if r46!=9: raise AssertionError(('T4 feedback rank',r46))
# The 49-dimensional physical complement is also independent under the diagonal K66 block.
colsW=[st['cols66'][i] for i in st['silent']]
rW=len(C.independent(colsW)[0])
if rW!=49: raise AssertionError(('silent K66 rank',rW))
# Sequential degree-six calibration supplies the natural first guess y_seq.
seq=C.solve_degree6_servo()
yseq=[seq['coeff'][i] for i in st['t4idx']]
# Diagnostic observer: pin to the midpoint of the sequential enclosure before nested inversions.
# Rank/structural reduction statements above do not depend on this midpoint calibration.
yobs=[C.arb(v.mid()) for v in yseq]
fb=H.feedback_map(st,yobs)
Ffield=C.combine(fb['F'],[st['V6b'][i] for i in st['t4idx']])
# Degree-four is exactly restored by construction for every y.
F4=C.vadd(st['N4'],C.vadd(C.combine(fb['c4'],st['K44']),fb['B4']))
if not C.norm2v(F4).contains(0): raise AssertionError(('Schur degree4 elimination',C.norm2v(F4)))
if not C.norm2v(fb['res6']).contains(0): raise AssertionError(('K66 inverse residual',C.norm2v(fb['res6'])))
out={
 'arb_precision_bits':BITS,'status':'PASS',
 'degree4_null_dimension':30,'degree6_null_dimension':58,
 'physical_feedback_sector':'degree-six toroidal l=4','feedback_dimension':9,
 'lower_silent_complement_dimension':49,
 'K46_rank_on_T4_feedback_sector':r46,
 'K46_silent_complement_structural_zero':all(C.norm2v(v).contains(0) for v in silent_back),
 'K66_rank_on_lower_silent_complement':rW,
 'sequential_guess_T4_feedback_mean_square':str(C.norm2v(C.combine(yobs,[st['V6b'][i] for i in st['t4idx']]))),
 'Schur_feedback_mismatch_at_sequential_midpoint_mean_square':str(C.norm2v(Ffield)),
 'degree4_residual_after_midpoint_observer_elimination':str(C.norm2v(F4)),
 'degree6_linear_inverse_residual_at_midpoint_observer':str(C.norm2v(fb['res6'])),
 'interpretation':'The apparent coupled 30+58 degree-four/degree-six maintenance problem collapses to a nine-dimensional physical feedback law because only the degree-six toroidal l=4 sector carries a degree-three harmonic Hodge companion and can therefore backreact onto degree four.  The remaining 49-dimensional sector P1+P3+P5+P7+T6 is exactly lower-silent and has rank 49 under the diagonal K66 response.  Since K44 and K66 are invertible, degree four and the 49 silent degree-six amplitudes can be eliminated uniquely for each nine-component T4 feedback y.  Coupled compatibility is therefore the intrinsic fixed-point equation y=Pi_T4 K66^{-1}[-R6(V4(y),y)], not an arbitrary 88-variable truncation.'}
print(json.dumps(out,indent=2,allow_nan=False))
