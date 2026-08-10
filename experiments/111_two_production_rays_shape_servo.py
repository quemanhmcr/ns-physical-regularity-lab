import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rt3=arb(3).sqrt()

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def outer(a,b): return tuple(tuple(a[i]*b[j] for j in range(3)) for i in range(3))
def madd(*Ms): return tuple(tuple(sum(M[i][j] for M in Ms) for j in range(3)) for i in range(3))
def mscale(c,M): return tuple(tuple(c*M[i][j] for j in range(3)) for i in range(3))
def contract(A,B): return sum(A[i][j]*B[i][j] for i in range(3) for j in range(3))
def mv(A,v): return tuple(sum(A[i][j]*v[j] for j in range(3)) for i in range(3))
def Mg(u,v,g): return madd(outer(u,v),outer(v,u),mscale(-g,outer(u,u)),mscale(-g,outer(v,v)))
def Aaxis(u): return tuple(tuple(u[i]*u[j]-(arb(1)/3 if i==j else arb(0)) for j in range(3)) for i in range(3))
def certify_one(x,label,tol='1e-30'):
 t=arb(tol)
 if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
n=(arb(0),arb(0),arb(1)); a=((arb(2)/3).sqrt(),arb(0),1/arb(3).sqrt()); b=(arb(0),-(arb(2)/3).sqrt(),-1/arb(3).sqrt())
Aa=Aaxis(a); Ab=Aaxis(b); An=Aaxis(n); Ma=Mg(a,n,dot(a,n)); Mb=Mg(b,n,dot(b,n)); Mc=Mg(a,b,dot(a,b))
rows=[]
for ss in ['1e-24','1','1e24']:
 s=arb(ss); la=arb(9)*s/4; lb=arb(9)*s/4
 Q2=madd(mscale(la,Aa),mscale(lb,Ab))
 ga=dot(a,mv(Q2,a)); gb=dot(b,mv(Q2,b)); certify_one(ga/s,('two-ray gain a',ss)); certify_one(gb/s,('two-ray gain b',ss))
 shape=[contract(M,Q2) for M in (Ma,Mb,Mc)]
 expected=[arb(2)*rt3*s/3,-arb(2)*rt3*s/3,-arb(4)*s/3]
 for i in range(3): certify_one(shape[i]/expected[i],('two-ray shape drift',ss,i))
 servo=mscale(-3*s,An)
 servo_gain_a=dot(a,mv(servo,a)); servo_gain_b=dot(b,mv(servo,b))
 if not (servo_gain_a.contains(0) and servo_gain_b.contains(0)): raise AssertionError(('servo should not change endpoint gains',ss,servo_gain_a,servo_gain_b))
 servo_shape=[contract(M,servo) for M in (Ma,Mb,Mc)]
 for i in range(3): certify_one(servo_shape[i]/(-expected[i]),('servo cancels shape drift',ss,i))
 locked=madd(Q2,servo)
 locked_shape=[contract(M,locked) for M in (Ma,Mb,Mc)]
 if not all(x.contains(0) for x in locked_shape): raise AssertionError(('servo-completed lock',ss,locked_shape))
 rows.append({'s':ss,'two_production_ray_coefficients':[str(la),str(lb)],'two_ray_shape_drift':[str(x) for x in shape],'servo_tensor_coefficient_minus3s':str(-3*s),'servo_endpoint_gains':[str(servo_gain_a),str(servo_gain_b)],'servo_shape_rates':[str(x) for x in servo_shape],'completed_shape_rates':[str(x) for x in locked_shape]})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'At the symmetric capacity pair, two axisymmetric endpoint-production rays can be chosen to give exactly s,s magnitude gains, but they necessarily drive nonzero Gram shape rates (2sqrt(3)s/3,-2sqrt(3)s/3,-4s/3). A third bridge-axis STF ray Q_servo=-3s A_n has exactly zero endpoint magnitude production and cancels all three shape drifts. Thus scalar mutual amplification does not close the tensor dynamics by itself: maintaining the productive chamber requires a pure shape-servo transaction channel.','rows':rows},indent=2,allow_nan=False))
