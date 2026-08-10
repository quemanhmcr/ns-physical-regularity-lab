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
  mm=arb(m); coeff=(mm-2)*(mm+3)
  for As in ['1e-24','1','1e24']:
    A=arb(As)
    for Ls in ['1e-12','1','1e12']:
      L=arb(Ls)
      for xs in ['1e-6','0.1','0.7']:
        x=arb(xs); r=x*L
        q=A*x**m
        qp=A*mm*r**(m-1)/L**m
        qpp=A*mm*(mm-1)*r**(m-2)/L**m
        L2=qpp+2*qp/r-6*q/(r*r)
        closed=coeff*A*x**(m-2)/(L*L)
        if m==2:
          if not L2.contains(0): raise AssertionError(('smooth l=2 zero mode lost',As,Ls,xs,L2))
        else:
          certify_one(L2/closed,('l=2 radial diffusion coefficient',m,As,Ls,xs))
        # Factorization L2 q = r^-4 d_r[r^6 d_r(q/r^2)].
        Kp=A*(mm-2)*r**(m-3)/L**m
        d_r_r6Kp=A*(mm-2)*(mm+3)*r**(m+2)/L**m
        fact=d_r_r6Kp/r**4
        if not (fact-L2).contains(0): raise AssertionError(('radial factorization',m,As,Ls,xs,fact,L2))
        rows.append({'m':m,'A':As,'L':Ls,'r_over_L':xs,'q':str(q),'L2_q':str(L2),'closed_coefficient':str(closed),'Q_over_r2_radial_derivative':str(Kp),'factorized_L2_q':str(fact)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'The sharp productive carrier is a toroidal l=2 field. Componentwise viscosity therefore acts on its STF radial transaction profile through L2 Q=Q_doubleprime+2Q_prime/r-6Q/r^2=r^-4 d_r[r^6 d_r(Q/r^2)]. The unique smooth radial viscosity-null profile is Q=r^2 C; this is exactly the tangent Hodge carrier profile. Local diffusion cannot be inferred from r alone when the productive carrier lies on this harmonic zero mode.','rows':rows},indent=2,allow_nan=False))
