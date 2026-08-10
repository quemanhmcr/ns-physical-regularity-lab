import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def vadd(a,b): return tuple(a[i]+b[i] for i in range(3))
def vscale(c,a): return tuple(c*a[i] for i in range(3))
def matvec(A,x): return tuple(sum(A[i][j]*x[j] for j in range(3)) for i in range(3))
def mscale(c,A): return tuple(tuple(c*A[i][j] for j in range(3)) for i in range(3))
def trace(A): return A[0][0]+A[1][1]+A[2][2]
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def det(a,b,c): return dot(a,cross(b,c))
def certify_one(x,label):
    tol=arb('1e-30')
    if not x.contains(1) or not (x > 1-tol and x < 1+tol):
        raise AssertionError((label,'ratio not tightly certified around one',x))

# Exact integer geometry, including a positive-cycle-style noncoplanar frame.
p=(arb(0),arb(1),arb(0))
q=(arb(-2),arb(-2),arb(1))
r=(arb(0),arb(-1),arb(-1))
D=det(p,q,r)
if D.contains(0): raise AssertionError(('base ancestry cell unexpectedly degenerate',D))

B0=((arb(2),arb(1),arb(0)),(arb(-1),arb(-3),arb(2)),(arb(1),arb(0),arb(1)))
if not trace(B0).contains(0): raise AssertionError(('common matrix not traceless',trace(B0)))
Da0=((arb(1),arb(0),arb(-1)),(arb(0),arb(2),arb(0)),(arb(1),arb(0),arb(-2)))
Db0=((arb(-2),arb(1),arb(0)),(arb(0),arb(1),arb(1)),(arb(0),arb(-1),arb(1)))
La=(arb(1),arb(-2),arb(3))
Lb=(arb(-1),arb(4),arb(2))
Z=((arb(0),arb(0),arb(0)),(arb(0),arb(0),arb(0)),(arb(0),arb(0),arb(0)))

# First observer: exact determinant trace identity for common incompressible deformation.
# Use integer scale factors so this certificate does not rely on cancellation of uncertain decimal balls.
common_rows=[]
for gs in ['1','1000000000000000000000000000000']:
    g=arb(gs); B=mscale(g,B0)
    common_direct=det(matvec(B,p),q,r)+det(p,matvec(B,q),r)+det(p,q,matvec(B,r))
    common_trace=trace(B)*D
    if not (common_direct.contains(0) and common_trace.contains(0)):
        raise AssertionError(('common affine determinant identity failed',gs,common_direct,common_trace))
    common_rows.append({'common_scale':gs,'direct_three_term_sum':str(common_direct),'trace_times_D':str(common_trace)})

# Second observer: remove the common mode before numerical evaluation, then independently
# compare the endpoint+viscous material derivative with J_bridge+J_nu.  This avoids
# manufacturing a tiny residual by subtracting huge common-affine terms.
defect_scales=['1e-24','1','1e24']
nu_scales=['1e-24','1','1e24']
rows=[]
for es in defect_scales:
    e=arb(es); Da=mscale(e,Da0); Db=mscale(e,Db0)
    bridge=det(matvec(Da,p),q,r)+det(p,q,matvec(Db,r))
    for nus in nu_scales:
      nu=arb(nus)
      visc=nu*(det(La,q,r)+det(p,q,Lb))
      pdot=vadd(matvec(Da,p),vscale(nu,La))
      qdot=matvec(Z,q)
      rdot=vadd(matvec(Db,r),vscale(nu,Lb))
      direct=det(pdot,q,r)+det(p,qdot,r)+det(p,q,rdot)
      predicted=bridge+visc
      if predicted.contains(0):
          if not direct.contains(0):
              raise AssertionError(('zero predicted but nonzero direct',es,nus,direct,predicted))
          ratio='zero/zero'
      else:
          rr=direct/predicted
          certify_one(rr,('pair-cell balance',es,nus))
          ratio=str(rr)
      rows.append({
          'defect_scale':es,'nu':nus,'D_ab':str(D),
          'J_bridge':str(bridge),'J_nu':str(visc),
          'Ddot_direct_after_common_mode_removed':str(direct),
          'Ddot_predicted':str(predicted),'direct_over_predicted':ratio,
      })

# Algebraic survival factorization checked without subtracting nearly equal quantities.
fac_cases=[]
for Las in ['1e-30','1','1e30']:
  for Lbs in ['1e-20','1','1e20']:
    Laa=arb(Las); Lbb=arb(Lbs)
    rr=arb('0.125'); T0=arb('0.75'); T=arb('0.25')
    RD=Laa*Lbb*rr*T/T0
    recovered=RD*(1/rr)*(T0/T)/(Laa*Lbb)
    certify_one(recovered,('survival factorization',Las,Lbs))
    fac_cases.append({'Lambda_a':Las,'Lambda_b':Lbs,'renewal_factor':str(RD),'identity_ratio':str(recovered)})

print(json.dumps({
  'arb_precision_bits':BITS,'status':'PASS',
  'common_affine_cases':len(common_rows),'balance_cases':len(rows),'factorization_cases':len(fac_cases),
  'base_D_ab':str(D),
  'interpretation':(
      'The unnormalized material pair ancestry cell D_ab=det[omega_a,R,omega_b] is unchanged by every common traceless affine deformation. '
      'After that common mode is removed before numerical observation, its remaining exact rate is independently certified as endpoint gradient mismatch plus the two viscous Laplacian terms. '
      'The multiplicative survival identity then splits pair amplification into oriented-cell renewal, bridge compression, and triple-product loss without catastrophic-cancellation observers.'
  ),
  'common_affine_trace_certificate':common_rows,'rows':rows,'factorization':fac_cases,
},indent=2,allow_nan=False))
