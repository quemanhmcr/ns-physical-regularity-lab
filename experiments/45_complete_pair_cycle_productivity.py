import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def det(a,b,c): return dot(a,cross(b,c))
def norm(a): return dot(a,a).sqrt()
def matvec(A,x): return tuple(sum(A[i][j]*x[j] for j in range(3)) for i in range(3))
def madd(A,B): return tuple(tuple(A[i][j]+B[i][j] for j in range(3)) for i in range(3))
def mscale(c,A): return tuple(tuple(c*A[i][j] for j in range(3)) for i in range(3))
def vadd(a,b): return tuple(a[i]+b[i] for i in range(3))
def vscale(c,a): return tuple(c*a[i] for i in range(3))
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x > 1-t and x < 1+t):
        raise AssertionError((label,'ratio not tightly certified around one',x))

def transpose(A): return tuple(tuple(A[j][i] for j in range(3)) for i in range(3))

# A positive mutual-stretching pair: D*L_a>0 and D*L_b<0.
p=(arb(2),arb(-1),arb(3))
R=(arb(1),arb(4),arb(-2))
q=(arb(-2),arb(0),arb(-3))
D=det(p,R,q); La=dot(p,R); Lb=dot(q,R)
if not (D*La>0 and D*Lb<0): raise AssertionError(('base pair not positive cycle',D,La,Lb))
ra=norm(p); rb=norm(q); rr=norm(R)
T=D/(ra*rb*rr); alpha=La/(ra*rr); beta=Lb/(rb*rr)
Kba=T*alpha; Kab=-T*beta; G=-alpha*beta; P=Kba*Kab
if not (Kba>0 and Kab>0 and G>0 and P>0): raise AssertionError('positive product gate failed')
certify_one(Kba/(D*La/(ra*ra*rb*rr*rr)),'edge b->a ancestry factorization')
certify_one(Kab/(-D*Lb/(ra*rb*rb*rr*rr)),'edge a->b ancestry factorization')
certify_one(P/(-D*D*La*Lb/(ra**3*rb**3*rr**4)),'cycle product ancestry factorization')
certify_one(P/(T*T*G),'P=T^2 G')

# Exact dynamic balance. Abar is incompressible common chord gradient; endpoint defects
# and viscous vectors are varied over extreme scales.  The common affine part must
# cancel from Ddot but act through Sbar in longitudinal overlap renewal.
Abar=((arb(2),arb(1),arb(0)),(arb(-1),arb(-1),arb(2)),(arb(3),arb(0),arb(-1)))
if not (Abar[0][0]+Abar[1][1]+Abar[2][2]).contains(0): raise AssertionError('Abar not tracefree')
At=transpose(Abar)
Sbar=tuple(tuple((Abar[i][j]+At[i][j])/2 for j in range(3)) for i in range(3))
Da0=((arb(1),arb(-2),arb(3)),(arb(0),arb(2),arb(-1)),(arb(4),arb(1),arb(-3)))
Db0=((arb(-2),arb(1),arb(0)),(arb(3),arb(-1),arb(2)),(arb(1),arb(4),arb(2)))
lap_p0=(arb(1),arb(-3),arb(2)); lap_q0=(arb(-2),arb(5),arb(1))
rows=[]
for ds in ['1e-24','1','1e24']:
  d=arb(ds)
  Da=mscale(d,Da0); Db=mscale(d,Db0)
  Aa=madd(Abar,Da); Ab=madd(Abar,Db)
  for ns in ['1e-24','1','1e24']:
    nu=arb(ns)
    lap_p=lap_p0; lap_q=lap_q0
    pdot=vadd(matvec(Aa,p),vscale(nu,lap_p))
    qdot=vadd(matvec(Ab,q),vscale(nu,lap_q))
    Rdot=matvec(Abar,R)
    Ddot=det(pdot,R,q)+det(p,Rdot,q)+det(p,R,qdot)
    Jbridge=det(matvec(Da,p),R,q)+det(p,R,matvec(Db,q))
    Jnu=nu*(det(lap_p,R,q)+det(p,R,lap_q))
    predD=Jbridge+Jnu
    certify_one(Ddot/predD,('D bridge/visc balance',ds,ns))

    Ladot=dot(pdot,R)+dot(p,Rdot)
    Lbdot=dot(qdot,R)+dot(q,Rdot)
    predLa=2*dot(p,matvec(Sbar,R))+dot(matvec(Da,p),R)+nu*dot(lap_p,R)
    predLb=2*dot(q,matvec(Sbar,R))+dot(matvec(Db,q),R)+nu*dot(lap_q,R)
    certify_one(Ladot/predLa,('La longitudinal balance',ds,ns))
    certify_one(Lbdot/predLb,('Lb longitudinal balance',ds,ns))

    sigma_a=dot(p,pdot)/(ra*ra); sigma_b=dot(q,qdot)/(rb*rb); sigma_R=dot(R,Rdot)/(rr*rr)
    lambda_D=Ddot/D; lambda_La=Ladot/La; lambda_Lb=Lbdot/Lb
    Tlog=lambda_D-sigma_a-sigma_b-sigma_R
    Glog=lambda_La+lambda_Lb-sigma_a-sigma_b-2*sigma_R
    Plog=2*lambda_D+lambda_La+lambda_Lb-3*sigma_a-3*sigma_b-4*sigma_R
    if not (Plog-(2*Tlog+Glog)).contains(0): raise AssertionError(('full cycle log decomposition',ds,ns,Plog,2*Tlog+Glog))
    rows.append({'defect_scale':ds,'nu':ns,'Ddot_over_D':str(lambda_D),'La_dot_over_La':str(lambda_La),'Lb_dot_over_Lb':str(lambda_Lb),'T_log_rate':str(Tlog),'G_log_rate':str(Glog),'P_log_rate':str(Plog),'D_balance_ratio':str(Ddot/predD),'La_balance_ratio':str(Ladot/predLa),'Lb_balance_ratio':str(Lbdot/predLb)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'base_pair':{'D':str(D),'L_a':str(La),'L_b':str(Lb),'T':str(T),'alpha':str(alpha),'beta':str(beta),'G_minus_alpha_beta':str(G),'K_b_to_a':str(Kba),'K_a_to_b':str(Kab),'P_edge_product':str(P)},
 'interpretation':(
   'The exact mutual-stretching product factorizes as P=K_ba K_ab=-T^2 alpha beta and equivalently through the unnormalized ancestry observables D, L_a=omega_a.R, L_b=omega_b.R. '
   'The dynamic calibration shows that common incompressible affine deformation cancels from D renewal but enters both longitudinal overlaps through the common strain Sbar. '
   'Consequently full cycle survival has two physically distinct gates and the exact logarithmic deficit decomposes as Delta_P=2 Delta_T+Delta_G.'
 ),
 'rows':rows,
},indent=2,allow_nan=False))
