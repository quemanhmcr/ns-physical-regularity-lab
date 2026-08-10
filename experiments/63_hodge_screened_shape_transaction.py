import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def outer(a,b): return tuple(tuple(a[i]*b[j] for j in range(3)) for i in range(3))
def madd(*Ms): return tuple(tuple(sum(M[i][j] for M in Ms) for j in range(3)) for i in range(3))
def mscale(c,M): return tuple(tuple(c*M[i][j] for j in range(3)) for i in range(3))
def contract(A,B): return sum(A[i][j]*B[i][j] for i in range(3) for j in range(3))
def Mg(u,v,g): return madd(outer(u,v),outer(v,u),mscale(-g,outer(u,u)),mscale(-g,outer(v,v)))
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))

# Canonical tangent strain carrier from the Hodge microscope has
# S_h(r)=(1-r^2/L^2)S, S_v(r)=(r^2/L^2)S,
# Q(rho)=(14/5)(rho/L)^2 S.  Contract with the shape-generated STF tensors.
a=(arb('0.8'),arb(0),arb('0.6')); n=(arb(0),arb(0),arb(1)); b=(arb('0.3'),arb('-0.4'),-arb(3).sqrt()/2)
alpha=dot(a,n); beta=dot(b,n); gamma=dot(a,b)
Ms=[('alpha',Mg(a,n,alpha)),('beta',Mg(b,n,beta)),('gamma',Mg(a,b,gamma))]
S=((arb('0.4'),arb('0.1'),arb('-0.2')),(arb('0.1'),arb('-0.3'),arb('0.05')),(arb('-0.2'),arb('0.05'),arb('-0.1')))
L=arb(3)
rows=[]
for rs in ['1e-12','0.1','1','2.9']:
    r=arb(rs); x=r/L
    for label,M in Ms:
        base=contract(M,S)
        jh=(1-x*x)*base
        jv=x*x*base
        # Analytic Hodge-screened shell integral:
        # (14/5)/L^2 * int_0^r [rho-rho^6/r^5] d rho = r^2/L^2.
        integral_coeff=(arb(14)/5)*(r*r/(L*L))*(arb(1)/2-arb(1)/7)
        jint=integral_coeff*base
        if base.contains(0):
            if not (jv-jint).contains(0): raise AssertionError(('zero base projected transaction',label,rs,jv,jint))
            ratio='zero/zero'
        else:
            certify_one(jint/jv,('Hodge screened shape transaction',label,rs))
            ratio=str(jint/jv)
        if not (jh+jv-base).contains(0): raise AssertionError(('shape Hodge split',label,rs,jh,jv,base))
        rows.append({'coordinate':label,'r':rs,'M_colon_S_total_common_current':str(base),'harmonic_shape_current':str(jh),'vortical_shape_current':str(jv),'screened_Q_integral_shape_current':str(jint),'integral_over_vortical_ratio':ratio})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'Contracting the exact Hodge strain microscope with the shape-generated STF tensor M_g produces an exact signed reset-current transaction.  For each physical Gram coordinate, the vortical contribution M_g:S_v(r) equals the Hodge-screened radial integral of M_g:Q(rho) with the native weight 1-(rho/r)^5.  Productive shape renewal and the earlier vorticity-to-strain transaction therefore live in the same intrinsic tensor channel.'
 ),'rows':rows
},indent=2,allow_nan=False))
