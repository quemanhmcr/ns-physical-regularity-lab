import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def certify_one(x,label,tol='1e-30'):
 t=arb(tol)
 if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
rows=[]
for As in ['1e-24','1','1e24']:
 A=arb(As); c4=arb('0.7'); c8=arb('-0.2')
 for nus in ['1e-24','1','1e24']:
  nu=arb(nus)
  for Ls in ['1e-12','1','1e12']:
   L=arb(Ls)
   for xs in ['1e-6','0.1','0.8']:
    x=arb(xs); r=x*L
    residual=A*(c4*x**4+c8*x**8)
    # K=transaction(curl j)=-nu L2 Q. Monomial L2 coefficients are 14 and 66.
    K4=-nu*A/L**2*arb(14)*c4
    K8=-nu*A/L**2*arb(66)*c8
    # -r^2/(5nu) int screen K ds/s, with K4*(s/L)^2 and K8*(s/L)^6.
    rec4=-(r*r/(5*nu))*K4*x*x*(arb(5)/(arb(2)*arb(7)))
    rec8=-(r*r/(5*nu))*K8*x**6*(arb(5)/(arb(6)*arb(11)))
    rec=rec4+rec8
    certify_one(rec/residual,('ancestry-current curvature reconstructs radial transaction defect',As,nus,Ls,xs))
    rows.append({'A':As,'nu':nus,'L':Ls,'r_over_L':xs,'radial_transaction_defect':str(residual),'curl_j_transaction_coefficients':[str(K4),str(K8)],'ancestry_curvature_reconstruction':str(rec),'ratio':str(rec/residual)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'Combining transaction(curl j)=-nu L2 Q with the l=2 Green identity gives Q-r^2 C=-(r^2/(5nu)) integral_0^r [1-(rho/r)^5] K_j(rho) d rho/rho, where K_j is the Hodge transaction tensor of curl j. The transaction profile that departs from the circulation-preserving smooth zero mode is therefore exactly a Hodge-screened accumulation of Kelvin ancestry-current curvature.','rows':rows},indent=2,allow_nan=False))
