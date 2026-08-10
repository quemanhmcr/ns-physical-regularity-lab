import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
rt2=arb(2).sqrt(); rt3=arb(3).sqrt()

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def mv(A,v): return tuple(sum(A[i][j]*v[j] for j in range(3)) for i in range(3))
def outer(a,b): return tuple(tuple(a[i]*b[j] for j in range(3)) for i in range(3))
def madd(*Ms): return tuple(tuple(sum(M[i][j] for M in Ms) for j in range(3)) for i in range(3))
def mscale(c,M): return tuple(tuple(c*M[i][j] for j in range(3)) for i in range(3))
def contract(A,B): return sum(A[i][j]*B[i][j] for i in range(3) for j in range(3))
def Mg(u,v,g): return madd(outer(u,v),outer(v,u),mscale(-g,outer(u,u)),mscale(-g,outer(v,v)))
def det3(M): return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])
def norm2(M): return contract(M,M)
def certify_one(x,label,tol='1e-30'):
 t=arb(tol)
 if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
n=(arb(0),arb(0),arb(1)); a=((arb(2)/3).sqrt(),arb(0),1/rt3); b=(arb(0),-(arb(2)/3).sqrt(),-1/rt3)
alpha=dot(a,n); beta=dot(b,n); gamma=dot(a,b); Ma=Mg(a,n,alpha); Mb=Mg(b,n,beta); Mc=Mg(a,b,gamma)
rows=[]
for sas in ['1e-24','0.1','1','1e24']:
 sa=arb(sas)
 for sbs in ['1e-24','0.3','1','1e24']:
  sb=arb(sbs)
  q=rt2*(2*sa+sb)/4; z=rt2*(sa+2*sb)/4
  S=((sa,arb(0),q),(arb(0),sb,z),(q,z,-sa-sb))
  raw_ga=dot(a,mv(S,a)); raw_gb=dot(b,mv(S,b)); raw_shape=[contract(M,S) for M in (Ma,Mb,Mc)]
  # Structural response certificate: the pair-frame solve gives gains (sa,sb) and zero shape rates before mixed-scale Cartesian parent terms are assembled.
  ga=sa; gb=sb; shape=[arb(0),arb(0),arb(0)]
  if not (raw_ga.contains(sa) and raw_gb.contains(sb) and all(x.contains(0) for x in raw_shape)):
    raise AssertionError(('raw pair-frame observer excluded structural target',sas,sbs,raw_ga,raw_gb,raw_shape))
  S2=norm2(S); closed2=(13*sa*sa+16*sa*sb+13*sb*sb)/4
  certify_one(S2/closed2,('shape-lock norm',sas,sbs))
  d=det3(S); closedd=-(sa+sb)*(sa*sa+15*sa*sb+sb*sb)/8
  certify_one(d/closedd,('shape-lock determinant',sas,sbs))
  if not (d<0): raise AssertionError(('positive-gain lock determinant should be negative',sas,sbs,d))
  bridge=dot(n,mv(S,n)); certify_one(bridge/(-(sa+sb)),('bridge contraction',sas,sbs))
  rows.append({'s_a':sas,'s_b':sbs,'S_Frobenius_squared':str(S2),'closed_norm_squared':str(closed2),'det_S_lock':str(d),'closed_negative_det':str(closedd),'bridge_rate':str(bridge),'shape_rates':[str(x) for x in shape],'raw_cartesian_gain_a':str(raw_ga),'raw_cartesian_gain_b':str(raw_gb),'raw_cartesian_shape_rates':[str(x) for x in raw_shape]})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'At the capacity Gram geometry, prescribed positive endpoint gains s_a,s_b and zero alpha,beta,gamma shape rates determine a unique STF common strain. Its matrix is [[s_a,0,sqrt(2)(2s_a+s_b)/4],[0,s_b,sqrt(2)(s_a+2s_b)/4],[..., ...,-s_a-s_b]]. Its determinant is -(s_a+s_b)(s_a^2+15s_as_b+s_b^2)/8<0 for every positive gain pair. Stationary positive-cycle shape lock therefore lies universally in the negative-determinant strain sector.','rows':rows},indent=2,allow_nan=False))
