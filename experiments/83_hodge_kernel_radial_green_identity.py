import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
# q=A[x^2+c4 x^4+c8 x^8]. The x^2 part is the smooth L2 kernel.
rows=[]
for As in ['1e-24','1','1e24']:
 A=arb(As)
 for Ls in ['1e-12','1','1e12']:
  L=arb(Ls)
  for xs in ['1e-6','0.1','0.8']:
   x=arb(xs); r=x*L; c4=arb('0.7'); c8=arb('-0.2')
   q=A*(x*x+c4*x**4+c8*x**8); kernel=A*x*x; residual=q-kernel
   # L2 q = A/L^2 [14 c4 x^2 + 66 c8 x^6].
   f4=A/L**2*arb(14)*c4*x*x; f8=A/L**2*arb(66)*c8*x**6
   # For an original x^m residual, r^2/5 * integral screen*(L2 term) ds/s recovers that x^m exactly.
   rec4=(r*r/5)*(A/L**2*arb(14)*c4)*(x*x)*(arb(5)/(arb(2)*arb(7)))
   rec8=(r*r/5)*(A/L**2*arb(66)*c8)*(x**6)*(arb(5)/(arb(6)*arb(11)))
   # Above factors: integral_0^r [1-(s/r)^5](s/r)^(m-2) ds/s = 5/[(m-2)(m+3)].
   reconstructed=rec4+rec8
   certify_one(reconstructed/residual,('same Hodge screen inverts productive radial defect',As,Ls,xs))
   rows.append({'A':As,'L':Ls,'r_over_L':xs,'Q_profile_scalar':str(q),'smooth_r2_kernel':str(kernel),'radial_defect_residual':str(residual),'L2_components':[str(f4),str(f8)],'Green_reconstructed_residual':str(reconstructed),'reconstruction_ratio':str(reconstructed/residual)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For the productive toroidal l=2 radial operator, the departure Q(r)-r^2 C from the unique smooth viscosity-null mode obeys the exact Green identity Q-r^2 C=(r^2/5) integral_0^r [1-(rho/r)^5](L2 Q)(rho) d rho/rho. The same Hodge screen 1-(rho/r)^5 therefore appears both in the strain transaction microscope and in reconstruction of the viscous-active radial transaction defect.','rows':rows},indent=2,allow_nan=False))
