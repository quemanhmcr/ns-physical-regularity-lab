import json, os, itertools
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def outer(a,b): return tuple(tuple(a[i]*b[j] for j in range(3)) for i in range(3))
def madd(*Ms): return tuple(tuple(sum(M[i][j] for M in Ms) for j in range(3)) for i in range(3))
def mscale(c,M): return tuple(tuple(c*M[i][j] for j in range(3)) for i in range(3))
def contract(A,B): return sum(A[i][j]*B[i][j] for i in range(3) for j in range(3))
def Aamp(u):
    return ((u[0]*u[0]-arb(1)/3,u[0]*u[1],u[0]*u[2]),(u[1]*u[0],u[1]*u[1]-arb(1)/3,u[1]*u[2]),(u[2]*u[0],u[2]*u[1],u[2]*u[2]-arb(1)/3))
def Mg(u,v,g): return madd(outer(u,v),outer(v,u),mscale(-g,outer(u,u)),mscale(-g,outer(v,v)))
def coord(M): return (M[0][0],M[1][1],M[0][1],M[0][2],M[1][2])
def parity(p):
    inv=sum(1 for i in range(len(p)) for j in range(i+1,len(p)) if p[i]>p[j])
    return -1 if inv%2 else 1
def det5(M):
    out=arb(0)
    for p in itertools.permutations(range(5)):
        term=arb(parity(p))
        for i in range(5): term*=M[i][p[i]]
        out+=term
    return out
def replace_col(A,j,y): return [[y[i] if k==j else A[i][k] for k in range(5)] for i in range(5)]
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))

a=(arb('0.8'),arb(0),arb('0.6')); n=(arb(0),arb(0),arb(1)); b=(arb('0.3'),arb('-0.4'),-arb(3).sqrt()/2)
alpha=dot(a,n); beta=dot(b,n); gamma=dot(a,b)
E=[Aamp(a),Aamp(b),Mg(a,n,alpha),Mg(b,n,beta),Mg(a,b,gamma)]
H=[[arb(2),arb(1),arb(0),arb(0),arb(0)],[arb(1),arb(2),arb(0),arb(0),arb(0)],[arb(0),arb(0),arb(2),arb(0),arb(0)],[arb(0),arb(0),arb(0),arb(2),arb(0)],[arb(0),arb(0),arb(0),arb(0),arb(2)]]
C=[[coord(E[j])[i] for j in range(5)] for i in range(5)]
Bmat=[[sum(C[k][i]*H[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
DB=det5(Bmat)
if DB.contains(0): raise AssertionError(('pair response matrix singular',DB))
rows=[]
for ss in ['1e-24','1','1e24']:
    q=arb(ss); svec=[arb('0.31')*q,-arb('0.27')*q,arb('0.19')*q,-arb('0.23')*q,arb('0.17')*q]
    Smat=((svec[0],svec[2],svec[3]),(svec[2],svec[1],svec[4]),(svec[3],svec[4],-svec[0]-svec[1]))
    responses=[contract(Ei,Smat) for Ei in E]
    recovered=[det5(replace_col(Bmat,j,responses))/DB for j in range(5)]
    for j in range(5): certify_one(recovered[j]/svec[j],('pair strain reconstruction',ss,j))
    rows.append({'strain_scale':ss,'responses_sigma_a_sigma_b_alpha_beta_gamma':[str(x) for x in responses],'original_STF_coordinates':[str(x) for x in svec],'recovered_STF_coordinates':[str(x) for x in recovered]})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'response_matrix_det':str(DB),
 'interpretation':(
  'For a noncoplanar productive pair, the five directly physical responses consisting of two magnitude-stretching rates and three Gram-shape rates reconstruct an arbitrary symmetric trace-free strain exactly.  The pair itself therefore supplies a complete observer frame without privileging Cartesian strain components.'
 ),'rows':rows
},indent=2,allow_nan=False))
