import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import degree6_hodge_servo_core as C
z=C.z;one=C.o
X=({(1,0,0):one},{(0,1,0):one},{(0,0,1):one})

def eval0(P): return P.get((0,0,0),z)
def gradmat(V): return tuple(tuple(eval0(C.pder(V[i],j)) for j in range(3)) for i in range(3))
def mm(A,B): return tuple(tuple(sum((A[i][k]*B[k][j] for k in range(3)),z) for j in range(3)) for i in range(3))
def madd(A,B): return tuple(tuple(A[i][j]+B[i][j] for j in range(3)) for i in range(3))
def mscale(c,A): return tuple(tuple(c*A[i][j] for j in range(3)) for i in range(3))
def msub(A,B): return madd(A,mscale(-1,B))
def mn2(A): return sum((A[i][j]*A[i][j] for i in range(3) for j in range(3)),z)
def tr(A): return sum((A[i][i] for i in range(3)),z)
def mpow(A,n):
 R=tuple(tuple(one if i==j else z for j in range(3)) for i in range(3))
 for _ in range(n):R=mm(R,A)
 return R

def matvec(A,V):
 out=[]
 for i in range(3):
  q={}
  for j in range(3):q=C.padd(q,C.pscale(A[i][j],V[j]))
  out.append(q)
 return tuple(out)
def pconst(c):return {(0,0,0):c}
def curlpot(Psi):return C.curl(Psi)

def null_u2(B): return C.vscale(-arb(1)/3,C.cross(X,matvec(B,X)))

Ss=[
 ((arb(2),z,z),(z,-one,z),(z,z,-one)),
 ((arb('0.4'),arb('0.2'),arb('-0.1')),(arb('0.2'),arb('-0.1'),arb('0.3')),(arb('-0.1'),arb('0.3'),arb('-0.3'))),
]
Bs=[
 ((arb(2),z,z),(z,-one,z),(z,z,-one)),
 ((one,arb('0.3'),z),(arb('0.3'),arb('-0.2'),arb('0.4')),(z,arb('0.4'),arb('-0.8'))),
]
# degree-5 vector potentials -> degree-4 divergence-free velocity -> cubic vorticity -> nonzero grad Delta omega at origin.
Psis=[
 ({(5,0,0):arb('0.07'),(2,3,0):arb('-0.11')},{(0,5,0):arb('0.05'),(1,1,3):arb('0.09')},{(0,0,5):arb('-0.04'),(3,0,2):arb('0.13')}),
 ({(4,1,0):arb('0.2'),(0,2,3):arb('-0.17')},{(1,4,0):arb('-0.08'),(3,0,2):arb('0.06')},{(2,1,2):arb('0.15'),(0,4,1):arb('-0.03')}),
]
rows=[]
for si,S in enumerate(Ss):
 if not tr(S).contains(0):raise AssertionError(('S trace',si,tr(S)))
 for bi,B in enumerate(Bs):
  if not tr(B).contains(0):raise AssertionError(('B trace',bi,tr(B)))
  for pi,Psi in enumerate(Psis):
   ulin=matvec(S,X);u2=null_u2(B);u4=curlpot(Psi);u=C.vadd(ulin,C.vadd(u2,u4));omega=C.curl(u)
   if any(not eval0(q).contains(0) for q in u):raise AssertionError(('u0',si,bi,pi))
   if any(not eval0(q).contains(0) for q in omega):raise AssertionError(('omega0',si,bi,pi))
   A0=gradmat(u);B0=gradmat(omega);lapw=tuple(C.plap(q) for q in omega);C3=gradmat(lapw)
   for nus in ('1e-20','1','1e20'):
    nu=arb(nus);F=C.vadd(C.bracket(omega,u),C.vscale(nu,lapw));Bdot=gradmat(F)
    pred=madd(msub(mm(A0,B0),mm(B0,A0)),mscale(nu,C3));err=mn2(msub(Bdot,pred))
    if not err.contains(0):raise AssertionError(('jet law',si,bi,pi,nus,err))
    I2dot=2*tr(mm(B0,Bdot));I2pred=2*nu*tr(mm(B0,C3))
    I3dot=3*tr(mm(mm(B0,B0),Bdot));I3pred=3*nu*tr(mm(mm(B0,B0),C3))
    if not (I2dot-I2pred).contains(0):raise AssertionError(('I2',si,bi,pi,nus,I2dot,I2pred))
    if not (I3dot-I3pred).contains(0):raise AssertionError(('I3',si,bi,pi,nus,I3dot,I3pred))
    rows.append({'S_case':si,'B_case':bi,'quartic_velocity_case':pi,'nu':nus,'jet_law_error_squared':str(err),'tr_B2':str(tr(mm(B0,B0))),'Dt_tr_B2':str(I2dot),'viscous_only_tr_B2_prediction':str(I2pred),'tr_B3':str(tr(mm(mm(B0,B0),B0))),'Dt_tr_B3':str(I3dot),'viscous_only_tr_B3_prediction':str(I3pred),'commutator_matrix_norm_squared':str(mn2(msub(mm(A0,B0),mm(B0,A0))))})
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'Differentiate the exact Navier-Stokes vorticity equation at a material point with u=0 and omega=0.  For B=grad omega and A=grad u, polynomial certification gives D_t B=[A,B]+nu grad Delta omega.  The potentially non-affine Euler term (grad A) omega vanishes exactly at a vorticity zero. '
  'Consequently every spectral trace tr(B^k) is changed only by viscosity: D_t tr(B^k)=k nu tr(B^(k-1) grad Delta omega); the Euler commutator drops out of the trace.  In Euler, the eigenvalues of the linear vorticity-gradient germ are material invariants even though the matrix can be strongly conjugated by strain. '
  'This is a local PDE ancestry law for the transaction-null linear catalyst.  Its strength/eigenvalues cannot be amplified by reusable Euler deformation at a vorticity-zero center; actual amplitude growth requires viscous third-derivative current or switching to a different material center/zero.  The next module converts the critical growth law into the forced sub-core viscous length.'),
 'rows':rows
},indent=2,allow_nan=False))
