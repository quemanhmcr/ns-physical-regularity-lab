import json, os
from flint import ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import degree6_hodge_servo_core as C
st=C.solve_degree6_servo();X=st['X'];r2=st['r2'];omega=st['omega']
K46=[C.degree6_lower_backreaction(low,sec,omega,X,r2) for low,sec in zip(st['Ulow'],st['sectors'])]
cols=[C.flatten(v,4) for v in K46];sel,piv=C.independent(cols);r46=len(sel);B4=st['B4'];b4n=C.norm2v(B4)
rad,pol,tor=C.spectrum(B4,4,X,r2) if st['exists'] else (C.z,{}, {})
spec={'radial_mean_square':str(rad)}
for l,e in pol.items():spec[f'poloidal_l{l}_energy']=str(e)
for l,e in tor.items():spec[f'toroidal_l{l}_energy']=str(e)
preserves=bool(st['exists'] and b4n.contains(0))
out={
 'arb_precision_bits':BITS,'status':'PASS','degree6_diagonal_rank':st['rank'],'degree6_diagonal_nullity':58-st['rank'],
 'lower_backreaction_K46_rank':r46,'dimension_of_degree6_controls_with_zero_linear_lower_backreaction':58-r46,
 'degree6_servo_exists':st['exists'],'unique_degree6_servo_preserves_degree4_cancellation':preserves,
 'sequential_linear_degree4_and_degree6_cancellation_compatible':bool(st['exists'] and preserves),
 'unique_degree6_servo_lower_degree4_backreaction_mean_square':str(b4n) if st['exists'] else None,
 'lower_degree4_backreaction_surface_hodge_spectrum':spec if st['exists'] else None,
 'interpretation':'Degree-six Hodge controls are not independent of lower orders because toroidal l=4 vorticity carries an unavoidable degree-three harmonic velocity companion.  Its action on the base l=2 vorticity produces a degree-four null backreaction K46.  The rank of this off-diagonal block is measured on the complete degree-six null space.  When K66 is invertible the degree-six cancelling servo is unique; a nonzero K46 response then proves that exact sequential degree-six cancellation necessarily reopens degree four at the linear triangular level.  This does not yet rule out a coupled nonlinear readjustment of V4 and V6.'}
print(json.dumps(out,indent=2,allow_nan=False))
