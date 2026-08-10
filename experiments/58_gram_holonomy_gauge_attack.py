import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def mm(A,B): return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def mv(A,v): return tuple(sum(A[i][j]*v[j] for j in range(3)) for i in range(3))
def Rz(th):
    c=th.cos(); s=th.sin()
    return ((c,-s,arb(0)),(s,c,arb(0)),(arb(0),arb(0),arb(1)))
def gram(a,n,b):
    alpha=dot(a,n); beta=dot(b,n); gamma=dot(a,b); T=dot(a,cross(n,b))
    F=1+2*alpha*beta*gamma-alpha*alpha-beta*beta-gamma*gamma
    return alpha,beta,gamma,T,F

def certify_same(x,y,label):
    if not (x-y).contains(0): raise AssertionError((label,x,y))

# A nontrivial closed Gram path is prescribed intrinsically by alpha(s), beta and
# a transverse angle psi(s).  Multiplying the entire lifted triad by an arbitrary
# common Q(s) leaves exactly the same Gram history but can change the final lift.
pi=arb.pi(); beta=-arb('0.5')
ss=['0','0.125','0.25','0.5','0.75','1']
angles=['0','1','1e12','1e24']
rows=[]
for angs in angles:
    A=arb(angs)
    for ss_ in ss:
        s=arb(ss_); phase=2*pi*s
        alpha=arb('0.4')+arb('0.1')*phase.sin()
        psi=arb('0.7')+arb('0.2')*phase.sin()
        sa=(1-alpha*alpha).sqrt(); sb=(1-beta*beta).sqrt()
        n=(arb(0),arb(0),arb(1))
        a=(sa,arb(0),alpha)
        b=(sb*psi.cos(),sb*psi.sin(),beta)
        g0=gram(a,n,b)
        Q=Rz(A*s)
        ar=mv(Q,a); nr=mv(Q,n); br=mv(Q,b)
        g1=gram(ar,nr,br)
        for k,name in enumerate(['alpha','beta','gamma','T','detG']): certify_same(g0[k],g1[k],('common lift gauge',angs,ss_,name))
        rows.append({'final_rotation_parameter':angs,'s':ss_,'alpha':str(g0[0]),'beta':str(g0[1]),'gamma':str(g0[2]),'T':str(g0[3]),'Gram_det':str(g0[4])})

# The Gram loop closes at s=0 and s=1, while the common lift can end at any angle.
closures=[]
for angs in angles:
    A=arb(angs)
    def state(s):
        phase=2*pi*s; alpha=arb('0.4')+arb('0.1')*phase.sin(); psi=arb('0.7')+arb('0.2')*phase.sin()
        sa=(1-alpha*alpha).sqrt(); sb=(1-beta*beta).sqrt(); n=(arb(0),arb(0),arb(1))
        a=(sa,arb(0),alpha); b=(sb*psi.cos(),sb*psi.sin(),beta)
        Q=Rz(A*s); return gram(mv(Q,a),mv(Q,n),mv(Q,b))
    g0=state(arb(0)); g1=state(arb(1))
    for k in range(5): certify_same(g0[k],g1[k],('closed Gram loop',angs,k))
    closures.append({'final_rotation_parameter':angs,'closed_Gram_alpha':str(g1[0]),'closed_Gram_beta':str(g1[1]),'closed_Gram_gamma':str(g1[2]),'closed_Gram_T':str(g1[3])})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'closure_cases':len(closures),
 'interpretation':(
  'A closed productive Gram history does not determine a lifted SO(3) holonomy.  The same nontrivial Gram loop can be multiplied by an arbitrary time-dependent common rotation and therefore end at any chosen common rotation without changing alpha, beta, gamma, T or det G.  Any rotation holonomy assigned using an external connection is gauge-dependent rather than an intrinsic productive cost.'
 ),'rows':rows,'closures':closures
},indent=2,allow_nan=False))
