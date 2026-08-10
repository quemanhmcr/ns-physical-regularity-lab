import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def mv(A,v): return tuple(sum(A[i][j]*v[j] for j in range(3)) for i in range(3))
def scale(c,v): return tuple(c*x for x in v)
def add(a,b): return tuple(a[i]+b[i] for i in range(3))
def sub(a,b): return tuple(a[i]-b[i] for i in range(3))
E=((arb('0.4'),arb('0.3'),-arb('0.2')),(arb('0.3'),-arb('0.1'),arb('0.25')),(-arb('0.2'),arb('0.25'),-arb('0.3')))
pts=[(arb(1),arb(0),arb(0)),(arb('0.3'),arb('0.4'),arb('0.75').sqrt()),(arb('-0.6'),arb('0.2'),arb('0.6').sqrt())]
rows=[]
for m in [2,4,8]:
 mm=arb(m)
 for As in ['1e-24','1','1e24']:
  A=arb(As)
  for nus in ['1e-24','1','1e24']:
   nu=arb(nus)
   for Ls in ['1e-6','1','1e6']:
    L=arb(Ls); r=L*arb('0.7'); q=A*(r/L)**m; qp=A*mm*r**(m-1)/L**m
    for p in pts:
     pn=dot(p,p).sqrt(); n=tuple(x/pn for x in p); En=mv(E,n); f=dot(n,En); PEn=sub(En,scale(f,n))
     jr=5*nu*q/r*f
     jt=scale((5*nu/3)*(qp+q/r),PEn)
     j=add(scale(jr,n),jt)
     if m==2:
      C=A/(L*L); x=scale(r,n); gradpsi=scale(5*nu*C,mv(E,x)); err=tuple(j[i]-gradpsi[i] for i in range(3))
      if not all(z.contains(0) for z in err): raise AssertionError(('zero-mode Kelvin current not exact gradient',As,nus,Ls,err))
    rows.append({'m':m,'A':As,'nu':nus,'L':Ls,'r_over_L':'0.7','q':str(q),'q_prime':str(qp),'structural_current_formula':'j_r=5 nu q/r (n.Qhat.n); j_t=(5nu/3)(q_prime+q/r) P_n E n','zero_mode_exact_gradient':m==2})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For the sharp minimum productive carrier omega=-(5/3)q(r)n cross E n, the Kelvin viscous current j=nu curl omega has exact radial/tangential components j_r=5nu(q/r)(n.E.n) and j_t=(5nu/3)(q_prime+q/r)P_n E n. On the unique smooth radial zero mode q=C r^2 this collapses to j=5nu C E x=grad[(5nu C/2)x.E.x], so j can be nonzero while every loop period and dj vanish: nonzero viscous current alone is not ancestry replacement.','rows':rows},indent=2,allow_nan=False))
