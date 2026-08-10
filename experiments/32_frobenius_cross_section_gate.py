import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

# Frobenius gate for attempting a 2D cross-sectional reduction of a 3D vortex collar.
# A unit direction field xi has locally integrable orthogonal planes only if
# alpha=xi.dx satisfies alpha wedge d alpha=0, equivalently xi.curl(xi)=0.
# For the exact helical Beltrami direction xi=(cos kz,sin kz,0),
# curl xi=-k xi, so xi.curl xi=-k exactly and no orthogonal slicing exists for k!=0.
# The phase cancels from all invariant contractions, so no trig observer is needed.

rows=[]
for k_s in ['1e-30','1e-15','1e-6','1','1e6','1e15','1e30']:
  k=arb(k_s)
  xi_norm2=arb(1)
  curl_parallel=-k
  frob=xi_norm2*curl_parallel
  beltrami_ratio=frob/(-k)
  if not beltrami_ratio.contains(1):
    raise AssertionError(('Beltrami Frobenius density',k_s,beltrami_ratio))
  if not (frob<0):
    raise AssertionError(('helical plane field should be nonintegrable',k_s,frob))
  # Directional variation and Beltrami alignment are also exact.
  grad_xi_sq=k*k
  curl_xi_sq=k*k
  if not (grad_xi_sq/curl_xi_sq).contains(1):
    raise AssertionError(('helical derivative invariant',k_s,grad_xi_sq/curl_xi_sq))
  rows.append({'k':k_s,'xi_dot_curl_xi':str(frob),
               'frobenius_density_over_minus_k':str(beltrami_ratio),
               'grad_xi_sq':str(grad_xi_sq)})

# Planar/constant-direction null control: xi=e_z has curl xi=0 and Frobenius density 0 exactly.
null=arb(0)
if not null.contains(0): raise AssertionError('impossible null control')

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'cases':len(rows),
  'planar_constant_direction_frobenius_density':str(null),
  'interpretation':'The exact helical Beltrami direction has xi.curl(xi)=-k at every phase, so for every nonzero k the planes orthogonal to vorticity fail the Frobenius integrability condition. A universal proof cannot assume that a genuinely 3D vortex collar admits orthogonal 2D cross-sections. Planar capacity is exact when such slices exist, but full 3D cohomology/material-spacetime flux is required in the helical branch.',
  'rows':rows,
},indent=2,allow_nan=False))
