import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def mv(A,v): return tuple(sum(A[i][j]*v[j] for j in range(3)) for i in range(3))
def scale(c,v): return tuple(c*x for x in v)
def sub(a,b): return tuple(a[i]-b[i] for i in range(3))
def norm(v): return dot(v,v).sqrt()
Q=((arb('0.4'),arb('0.3'),-arb('0.2')),(arb('0.3'),-arb('0.1'),arb('0.25')),(-arb('0.2'),arb('0.25'),-arb('0.3')))
pts=[(arb(1),arb(0),arb(0)),(arb(0),arb(1),arb(0)),(arb('0.3'),arb('0.4'),arb('0.75').sqrt()),(arb('-0.6'),arb('0.2'),arb('0.6').sqrt())]
rows=[]
for p in pts:
    nn=norm(p); n=tuple(x/nn for x in p); Qn=mv(Q,n); qnn=dot(n,Qn)
    grad=scale(2,sub(Qn,scale(qnn,n))); omega=scale(-arb(5)/3,cross(n,Qn))
    raw_n=dot(n,omega); raw_g=dot(grad,omega)
    if not raw_n.contains(0) or not raw_g.contains(0): raise AssertionError(('raw structural cancellation excluded zero',raw_n,raw_g))
    # Structural certificate: n.(n cross Qn)=0 and [Qn-(n.Qn)n].(n cross Qn)=0 before parent terms are formed.
    structural_n=arb(0); structural_g=arb(0)
    rows.append({'n':[str(x) for x in n],'f_nQn':str(qnn),'gradS_f_norm':str(norm(grad)),'omega_ext_norm':str(norm(omega)),'structural_n_dot_omega':str(structural_n),'structural_grad_f_dot_omega':str(structural_g),'raw_n_dot_omega':str(raw_n),'raw_grad_f_dot_omega':str(raw_g)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'structural_surface_divergence':'zero: omega_ext=-(5/6)n cross grad_S f is a rotated surface gradient','interpretation':'The sharp equality field is tangent to the sphere and to level curves of f(n)=n.Q.n. The minimum transaction carrier therefore has intrinsic quadratic-level-set organization rather than an imposed filament closure.','rows':rows},indent=2,allow_nan=False))
