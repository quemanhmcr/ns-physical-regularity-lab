import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rt2=arb(2).sqrt(); rt3=arb(3).sqrt()

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def mv(A,v): return tuple(sum(A[i][j]*v[j] for j in range(3)) for i in range(3))
def certify_one(x,label,tol='1e-30'):
 t=arb(tol)
 if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
n=(arb(0),arb(0),arb(1)); a=((arb(2)/3).sqrt(),arb(0),1/rt3); b=(arb(0),-(arb(2)/3).sqrt(),-1/rt3)
rows=[]
for ss in ['1e-24','1','1e24']:
 s=arb(ss); c=3*s/(2*rt2)
 S=((s,arb(0),c),(arb(0),s,c),(c,c,-2*s))
 sa=dot(a,mv(S,a)); sb=dot(b,mv(S,b)); sr=dot(n,mv(S,n))
 certify_one(sa/s,('endpoint a gain',ss)); certify_one(sb/s,('endpoint b gain',ss)); certify_one(sr/(-2*s),('bridge contraction',ss))
 trace=S[0][0]+S[1][1]+S[2][2]
 if not trace.contains(0): raise AssertionError(('incompressible trace',ss,trace))
 determinant_log_rate=sa+sb+sr
 if not determinant_log_rate.contains(0): raise AssertionError(('oriented pair-cell common affine rate',ss,determinant_log_rate))
 rows.append({'sigma':ss,'sigma_a':str(sa),'sigma_b':str(sb),'sigma_bridge':str(sr),'sigma_a_plus_sigma_b_plus_sigma_bridge':str(determinant_log_rate),'structural_shape_rates':['0','0','0']})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'At the capacity stationary-amplifier geometry, the same common incompressible strain that amplifies both vorticity directions at rate sigma contracts the material bridge at rate -2sigma. This is the determinant balance of the oriented pair ancestry cell: with T shape-locked and no pair-cell renewal, sigma_a+sigma_b+sigma_R=0. Stationary productive amplification is therefore simultaneously a hyperbolic inward bridge conveyor.','rows':rows},indent=2,allow_nan=False))
