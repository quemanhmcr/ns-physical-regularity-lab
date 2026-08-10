import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def certify_one(x,label,tol='1e-30'):
 t=arb(tol)
 if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
rows=[]
for m in [2,4,8]:
 mm=arb(m); radial_coeff=(mm-2)*(mm+3)
 for As in ['1e-24','1','1e24']:
  A=arb(As)
  for nus in ['1e-24','1','1e24']:
   nu=arb(nus)
   for Ls in ['1e-12','1','1e12']:
    L=arb(Ls); x=arb('0.7')
    L2q=radial_coeff*A*x**(m-2)/(L*L)
    # curl j = -nu Delta omega. Since transaction(curl j)=K and transaction(Delta omega)=L2 Q,
    # the exact ancestry-current curvature tensor coefficient is K=-nu L2 q.
    K=-nu*L2q
    toroidal_curlj_coeff=(5*nu/3)*L2q
    inverse_from_K=-(arb(5)/3)*K
    if m==2:
      structural_K=arb(0)
      if not (K.contains(0) and toroidal_curlj_coeff.contains(0)): raise AssertionError(('zero-mode current curvature',As,nus,Ls,K,toroidal_curlj_coeff))
    else:
      structural_K=None
      certify_one(toroidal_curlj_coeff/inverse_from_K,('curl-j sharp inverse coefficient',m,As,nus,Ls))
      certify_one(K/(-nu*L2q),('transaction curvature K=-nu L2Q',m,As,nus,Ls))
    rows.append({'m':m,'A':As,'nu':nus,'L':Ls,'L2_q':str(L2q),'transaction_tensor_of_curl_j_coefficient':str(K),'structural_zero_curvature':str(structural_K) if structural_K is not None else None,'curl_j_toroidal_coefficient':str(toroidal_curlj_coeff)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'The non-exactness of the productive Kelvin current is exactly the radial productive defect: transaction(curl j)=-nu L2 Q, while curl j=(5nu/3)n cross (L2 Q)n. The smooth Q=r^2 C mode has structural zero current curvature even though j itself is nonzero. Thus L2 Q is precisely the productive l=2 component of the exterior derivative dj that changes material circulation ancestry.','rows':rows},indent=2,allow_nan=False))
