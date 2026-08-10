import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def norm(a): return dot(a,a).sqrt()
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x > 1-t and x < 1+t): raise AssertionError((label,x))
def certify_near_one(x,label,tol='1e-20'):
    t=arb(tol)
    if not (x > 1-t and x < 1+t): raise AssertionError((label,x))

theta_sets=[('-2','-0.3','1.1'),('-1','0.2','2'),('0.1','1','3')]
Ss=['1e-12','1','1e12','1e60']
eps=arb('3')
rows=[]; asym=[]
for ths in theta_sets:
  th=[arb(s) for s in ths]
  y=[t.cosh() for t in th]; z=[t.sinh() for t in th]
  c=[y[i]*z[i] for i in range(3)]
  C=[(2*th[i]).cosh() for i in range(3)]
  for Ss_ in Ss:
    S=arb(Ss_)
    X=[(S*c[i],y[i],z[i]) for i in range(3)]
    w=[(arb(0),eps*y[i],-eps*z[i]) for i in range(3)]
    # Oriented triangle closure remains exact while every unordered pair is productive.
    R12=tuple(X[1][k]-X[0][k] for k in range(3)); R23=tuple(X[2][k]-X[1][k] for k in range(3)); R31=tuple(X[0][k]-X[2][k] for k in range(3))
    for k in range(3):
      if not (R12[k]+R23[k]+R31[k]).contains(0): raise AssertionError(('vector closure',ths,Ss_,k))
    for i,j in [(0,1),(0,2),(1,2)]:
      R=tuple(X[j][k]-X[i][k] for k in range(3)); r=norm(R)
      rhoi=norm(w[i]); rhoj=norm(w[j]); n=tuple(v/r for v in R)
      xi=tuple(v/rhoi for v in w[i]); xj=tuple(v/rhoj for v in w[j])
      D=dot(w[i],cross(R,w[j])); Li=dot(w[i],R); Lj=dot(w[j],R)
      T=dot(xi,cross(n,xj)); alpha=dot(xi,n); beta=dot(xj,n)
      Kji=T*alpha; Kij=-T*beta; P=Kji*Kij
      delta=th[j]-th[i]; h=delta.cosh()-1; q=S*(c[j]-c[i]); r2=dot(R,R)
      Dcl=eps*eps*q*delta.sinh(); Licl=eps*h; Ljcl=-eps*h
      Tcl=q*delta.sinh()/((C[i]*C[j]).sqrt()*r)
      Gcl=h*h/((C[i]*C[j]).sqrt()*r2)
      CC=C[i]*C[j]; CC32=CC*CC.sqrt()
      Pcl=q*q*delta.sinh()**2*h*h/(CC32*r2*r2)
      certify_one(D/Dcl,('D closed',ths,Ss_,i,j)); certify_one(Li/Licl,('Li hyperbola',ths,Ss_,i,j)); certify_one(Lj/Ljcl,('Lj hyperbola',ths,Ss_,i,j))
      certify_one(T/Tcl,('T closed',ths,Ss_,i,j)); certify_one(P/Pcl,('P closed',ths,Ss_,i,j))
      if not (T>0 and alpha>0 and beta<0 and Kji>0 and Kij>0 and P>0):
        raise AssertionError(('positive pair cycle lost',ths,Ss_,i,j,T,alpha,beta,Kji,Kij,P))
      rows.append({'theta':ths,'S':Ss_,'pair':[i,j],'q':str(q),'D':str(D),'L_i':str(Li),'L_j':str(Lj),'T':str(T),'G':str(-alpha*beta),'K_j_to_i':str(Kji),'K_i_to_j':str(Kij),'P':str(P)})
      if Ss_=='1e60':
        Tinf=delta.sinh()/(C[i]*C[j]).sqrt()
        Pasym=delta.sinh()**2*h*h/(CC32*q*q)
        certify_near_one(T/Tinf,('T nonzero asymptote',ths,i,j))
        certify_near_one(P/Pasym,('P longitudinal starvation',ths,i,j))
        asym.append({'theta':ths,'pair':[i,j],'T_over_limit':str(T/Tinf),'P_over_qminus2_asymptote':str(P/Pasym),'T_limit':str(Tinf)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'asymptotic_cases':len(asym),
 'interpretation':(
  'Three material nodes chosen on y^2-z^2=1 form a closed triangle for which every unordered pair is a positive mutual-stretching two-cycle at each stroboscopic return. '
  'The common memory clock S drives every axial bridge like S while each transverse triple-product T approaches a nonzero limit. '
  'Nevertheless every full pair-cycle product decays like S^-2 because the opposite longitudinal accesses decay like 1/S. Closure does not reset productive memory.'
 ),'rows':rows,'asymptotic':asym
},indent=2,allow_nan=False))
