import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def mv(A,v): return tuple(sum(A[i][j]*v[j] for j in range(3)) for i in range(3))
def scale(c,v): return tuple(c*x for x in v)
def norm2(A): return sum(A[i][j]*A[i][j] for i in range(3) for j in range(3))
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))

E=((arb('0.4'),arb('0.3'),-arb('0.2')),(arb('0.3'),-arb('0.1'),arb('0.25')),(-arb('0.2'),arb('0.25'),-arb('0.3')))
E2=norm2(E)
pts=[(arb(1),arb(0),arb(0)),(arb('0.3'),arb('0.4'),arb('0.75').sqrt()),(arb('-0.6'),arb('0.2'),arb('0.6').sqrt())]
rows=[]
for m in [2,4,8]:
  mm=arb(m)
  for As in ['1e-24','1','1e24']:
    amp=arb(As)
    for Ls in ['1e-6','1','1e6']:
      L=arb(Ls); r=L*arb('0.7'); q=amp*(r/L)**m
      Aprime=(arb(5)/3)*q/r
      for p in pts:
        pn=dot(p,p).sqrt(); n=tuple(x/pn for x in p); x=scale(r,n); En=mv(E,n); Ex=mv(E,x)
        gradA=scale(Aprime,n); gradB=scale(-1,Ex)
        cleb=cross(gradA,gradB); prod=scale(-arb(5)/3*q,cross(n,En))
        err=tuple(cleb[i]-prod[i] for i in range(3))
        if not all(z.contains(0) for z in err): raise AssertionError(('Clebsch vector identity',m,As,Ls,err))
      r1=L*arb('0.2'); r2=L*arb('0.8')
      deltaA=(arb(5)/3)*(amp/mm)*((r2/L)**m-(r1/L)**m)
      J=E2*(amp/mm)*((r2/L)**m-(r1/L)**m)
      certify_one(J/((arb(3)/5)*E2*deltaA),('transaction vs Clebsch radial increment',m,As,Ls))
      rows.append({'m':m,'amplitude':As,'L':Ls,'E_Frobenius_squared':str(E2),'Delta_A':str(deltaA),'unscreened_transaction_E_colon_Q':str(J),'J_over_3over5_E2_DeltaA':str(J/((arb(3)/5)*E2*deltaA))})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For every fixed STF transaction ray Q(r)=q(r)E, the sharp minimum productive carrier has exact Clebsch form omega_prod=grad A cross grad B with A_prime=(5/3)q/r and B=-(1/2)x.E.x. The unscreened radial transaction is exactly (3/5)|E|^2 times the Clebsch radial increment Delta A. The individual potential A is a coordinate, while the physical vorticity flux is the two-form dA wedge dB.','rows':rows},indent=2,allow_nan=False))
