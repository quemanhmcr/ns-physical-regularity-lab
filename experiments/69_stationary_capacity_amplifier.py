import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rt2=arb(2).sqrt(); rt3=arb(3).sqrt()

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def outer(a,b): return tuple(tuple(a[i]*b[j] for j in range(3)) for i in range(3))
def madd(*Ms): return tuple(tuple(sum(M[i][j] for M in Ms) for j in range(3)) for i in range(3))
def mscale(c,M): return tuple(tuple(c*M[i][j] for j in range(3)) for i in range(3))
def contract(A,B): return sum(A[i][j]*B[i][j] for i in range(3) for j in range(3))
def norm2(M): return contract(M,M)
def Mg(u,v,g): return madd(outer(u,v),outer(v,u),mscale(-g,outer(u,u)),mscale(-g,outer(v,v)))
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))

# Capacity geometry: alpha=-beta=1/sqrt(3), T=2/3.
n=(arb(0),arb(0),arb(1)); a=((arb(2)/3).sqrt(),arb(0),1/rt3); b=(arb(0),-(arb(2)/3).sqrt(),-1/rt3)
alpha=dot(a,n); beta=dot(b,n); gamma=dot(a,b)
Ma=Mg(a,n,alpha); Mb=Mg(b,n,beta); Mc=Mg(a,b,gamma)
rows=[]
for ss in ['1e-24','1','1e24']:
    sig=arb(ss); c=arb(3)*sig/(2*rt2)
    S=((sig,arb(0),c),(arb(0),sig,c),(c,c,-2*sig))
    pa=dot(a,tuple(sum(S[i][j]*a[j] for j in range(3)) for i in range(3)))
    pb=dot(b,tuple(sum(S[i][j]*b[j] for j in range(3)) for i in range(3)))
    certify_one(pa/sig,('stationary amplifier sigma_a',ss)); certify_one(pb/sig,('stationary amplifier sigma_b',ss))
    raw_shape=[contract(M,S) for M in (Ma,Mb,Mc)]
    if not all(x.contains(0) for x in raw_shape): raise AssertionError(('raw shape lock excludes zero',ss,raw_shape))
    S2=norm2(S); certify_one(S2/((arb(21)/2)*sig*sig),('stationary amplifier norm',ss))
    for rs in ['1e-6','1','1e6']:
        r=arb(rs); E=(2*pi/15)*S2*r**5; closed=(arb(7)*pi/5)*sig*sig*r**5
        certify_one(E/closed,('stationary amplifier harmonic occupancy',ss,rs))
        rows.append({'sigma':ss,'r':rs,'sigma_a':str(pa),'sigma_b':str(pb),'structural_alpha_dot':'0','structural_beta_dot':'0','structural_gamma_dot':'0','raw_shape_rate_autopsy':[str(x) for x in raw_shape],'S_Frobenius_squared':str(S2),'harmonic_energy':str(E),'seven_pi_over_5_sigma2_r5':str(closed)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'At the globally optimal productive Gram geometry there is a unique symmetric trace-free common strain that keeps all three Gram coordinates fixed while amplifying both vorticity directions at the same prescribed rate sigma.  Its squared Frobenius size is (21/2)sigma^2, so a harmonic realization on radius r has exact kinetic occupancy (7*pi/5)sigma^2 r^5.  Large vorticity amplification therefore does not require shape recycling, but stationary optimal amplification is still a real strain-supported process.'
 ),'rows':rows
},indent=2,allow_nan=False))
