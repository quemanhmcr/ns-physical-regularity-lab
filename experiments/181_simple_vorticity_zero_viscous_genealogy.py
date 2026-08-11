import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import degree6_hodge_servo_core as C
z=C.z; one=C.o
X=({(1,0,0):one},{(0,1,0):one},{(0,0,1):one})

def pconst(a): return {(0,0,0):a}
def eval0(P): return P.get((0,0,0),z)
def gradmat(V): return tuple(tuple(eval0(C.pder(V[i],j)) for j in range(3)) for i in range(3))
def mm(A,B): return tuple(tuple(sum((A[i][k]*B[k][j] for k in range(3)),z) for j in range(3)) for i in range(3))
def madd(A,B): return tuple(tuple(A[i][j]+B[i][j] for j in range(3)) for i in range(3))
def mscale(c,A): return tuple(tuple(c*A[i][j] for j in range(3)) for i in range(3))
def msub(A,B): return madd(A,mscale(-1,B))
def mn2(A): return sum((A[i][j]*A[i][j] for i in range(3) for j in range(3)),z)
def tr(A): return sum((A[i][i] for i in range(3)),z)
def matvec(A,v): return tuple(sum((A[i][j]*v[j] for j in range(3)),z) for i in range(3))
def det3(A):
 a,b,c=A[0];d,e,f=A[1];g,h,i=A[2]
 return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
def inv3(A):
 a,b,c=A[0];d,e,f=A[1];g,h,i=A[2];D=det3(A)
 if D.contains(0): raise AssertionError(('singular B',D))
 return (( (e*i-f*h)/D,(c*h-b*i)/D,(b*f-c*e)/D),
         ( (f*g-d*i)/D,(a*i-c*g)/D,(c*d-a*f)/D),
         ( (d*h-e*g)/D,(b*g-a*h)/D,(a*e-b*d)/D))
def poly_matvec(A,V):
 out=[]
 for i in range(3):
  q={}
  for j in range(3): q=C.padd(q,C.pscale(A[i][j],V[j]))
  out.append(q)
 return tuple(out)
def null_u2(B): return C.vscale(-arb(1)/3,C.cross(X,poly_matvec(B,X)))
def quad_vort_velocity(qx,qy,qz):
 x,y,zz=X
 return (C.padd(C.pscale(qy/3,C.pmul(C.pmul(zz,zz),zz)),C.pscale(-qz/3,C.pmul(C.pmul(y,y),y))),
         C.padd(C.pscale(-qx/3,C.pmul(C.pmul(zz,zz),zz)),C.pscale(qz/3,C.pmul(C.pmul(x,x),x))),
         C.padd(C.pscale(qx/3,C.pmul(C.pmul(y,y),y)),C.pscale(-qy/3,C.pmul(C.pmul(x,x),x))))

def quartic_velocity(Psi): return C.curl(Psi)

def vconst(v): return tuple(pconst(a) for a in v)
def v0(V): return tuple(eval0(q) for q in V)
def vn2(v): return sum((a*a for a in v),z)
def trace_power_deriv(B,Bdot,k):
 P=tuple(tuple(one if i==j else z for j in range(3)) for i in range(3))
 for _ in range(k-1): P=mm(P,B)
 return arb(k)*tr(mm(P,Bdot))

Ss=[
 ((arb(2),z,z),(z,-one,z),(z,z,-one)),
 ((arb('0.4'),arb('0.2'),arb('-0.1')),(arb('0.2'),arb('-0.1'),arb('0.3')),(arb('-0.1'),arb('0.3'),arb('-0.3'))),
]
Bs=[
 ((arb(2),z,z),(z,-one,z),(z,z,-one)),
 ((one,arb('0.3'),z),(arb('0.3'),arb('-0.2'),arb('0.4')),(z,arb('0.4'),arb('-0.8'))),
]
qvecs=[(arb('0.07'),arb('-0.11'),arb('0.13')),(arb('-0.2'),arb('0.05'),arb('0.09'))]
Psis=[
 ({(5,0,0):arb('0.07'),(2,3,0):arb('-0.11')},{(0,5,0):arb('0.05'),(1,1,3):arb('0.09')},{(0,0,5):arb('-0.04'),(3,0,2):arb('0.13')}),
 ({(4,1,0):arb('0.2'),(0,2,3):arb('-0.17')},{(1,4,0):arb('-0.08'),(3,0,2):arb('0.06')},{(2,1,2):arb('0.15'),(0,4,1):arb('-0.03')}),
]
rows=[]
for si,S in enumerate(Ss):
 for bi,Btarget in enumerate(Bs):
  for qi,qv in enumerate(qvecs):
   for pi,Psi in enumerate(Psis):
    ulin=poly_matvec(S,X);u2=null_u2(Btarget);u3=quad_vort_velocity(*qv);u4=quartic_velocity(Psi)
    u=C.vadd(ulin,C.vadd(u2,C.vadd(u3,u4)));omega=C.curl(u)
    if vn2(v0(u)).contains(0) is False or vn2(v0(omega)).contains(0) is False: raise AssertionError(('origin',si,bi,qi,pi))
    A0=gradmat(u);B0=gradmat(omega);Binv=inv3(B0);lapw=tuple(C.plap(q) for q in omega);L0=v0(lapw);C3=gradmat(lapw)
    for nus in ('1e-20','1','1e20'):
     nu=arb(nus);w=tuple(-nu*a for a in matvec(Binv,L0));wv=vconst(w)
     F=C.vadd(C.bracket(omega,u),C.vscale(nu,lapw))
     zero_res=tuple(v0(F)[i]+matvec(B0,w)[i] for i in range(3))
     if not vn2(zero_res).contains(0): raise AssertionError(('zero drift law',si,bi,qi,pi,nus,zero_res))
     Hw=gradmat(C.directional(wv,omega))
     Bzdot=gradmat(C.vadd(F,C.directional(wv,omega)))
     pred=madd(msub(mm(A0,B0),mm(B0,A0)),madd(mscale(nu,C3),Hw))
     err=mn2(msub(Bzdot,pred))
     if not err.contains(0): raise AssertionError(('zero worldline B law',si,bi,qi,pi,nus,err))
     visc_eff=madd(mscale(nu,C3),Hw)
     for k in (2,3):
      lhs=trace_power_deriv(B0,Bzdot,k);rhs=trace_power_deriv(B0,visc_eff,k)
      if not (lhs-rhs).contains(0): raise AssertionError(('spectral zero genealogy',si,bi,qi,pi,nus,k,lhs,rhs))
     rows.append({'S_case':si,'B_case':bi,'quadratic_vorticity_case':qi,'cubic_vorticity_case':pi,'nu':nus,'det_grad_omega_at_zero':str(det3(B0)),'Delta_omega_at_zero':[str(x) for x in L0],'zero_relative_velocity_Vminus_u':[str(x) for x in w],'zero_condition_residual_square':str(vn2(zero_res)),'zero_worldline_B_law_error_squared':str(err),'Dt_zero_tr_B2':str(trace_power_deriv(B0,Bzdot,2)),'viscous_only_zero_tr_B2_prediction':str(trace_power_deriv(B0,visc_eff,2)),'Dt_zero_tr_B3':str(trace_power_deriv(B0,Bzdot,3)),'viscous_only_zero_tr_B3_prediction':str(trace_power_deriv(B0,visc_eff,3))})
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'Track a simple vorticity zero z(t), omega(z(t),t)=0, with invertible B=grad omega.  Differentiating the exact Navier-Stokes zero condition gives the zero-set drift V_z-u=-nu B^{-1} Delta omega.  Thus a simple vorticity zero is exactly material in Euler and can move relative to the fluid only through viscosity. '
  'Following that zero worldline gives dB/dt=[A,B]+nu grad Delta omega +(V_z-u).grad B.  Substituting the drift law shows that every change of tr(B^k) is viscosity-mediated; the Euler commutator again drops out.  Polynomial fields with nonzero Delta omega, nonzero grad Delta omega, noncommuting A and B, and two invertible B geometries certify both identities. '
  'This closes the fresh-simple-zero switching loophole at the local PDE level: changing which simple zero is followed does not create a new Euler source of catalyst spectral strength.  The remaining escape is loss of transversality det B->0, where zero branches can cease to be isolated and the inverse-B genealogy law becomes singular.'),
 'rows':rows
},indent=2,allow_nan=False))
