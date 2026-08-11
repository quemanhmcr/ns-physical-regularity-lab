import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import degree6_hodge_servo_core as C
z=C.z; one=C.o
X=({(1,0,0):one},{(0,1,0):one},{(0,0,1):one})

def matvec(B,V):
    out=[]
    for i in range(3):
        q={}
        for j in range(3): q=C.padd(q,C.pscale(B[i][j],V[j]))
        out.append(q)
    return tuple(out)

def mprod(A,B): return tuple(tuple(sum((A[i][k]*B[k][j] for k in range(3)),z) for j in range(3)) for i in range(3))
def tr(A): return sum((A[i][i] for i in range(3)),z)
def mtf(A):
    t=tr(A)/3
    return tuple(tuple(A[i][j]-(t if i==j else z) for j in range(3)) for i in range(3))
def mscale(c,A): return tuple(tuple(c*A[i][j] for j in range(3)) for i in range(3))
def msub(A,B): return tuple(tuple(A[i][j]-B[i][j] for j in range(3)) for i in range(3))
def mnorm2(A): return sum((A[i][j]*A[i][j] for i in range(3) for j in range(3)),z)
def vnorm2(V): return C.savg(C.vdot(V,V))
def lapv(V): return tuple(C.plap(v) for v in V)

def transaction_Q(V):
    nx=C.cross(X,V); Q=[]
    for i in range(3):
        row=[]
        for j in range(3):
            # Q=(3/2) sphere_avg[n_i (n x omega)_j + (n x omega)_i n_j]
            val=arb(3)/2*C.savg(C.padd(C.pmul(X[i],nx[j]),C.pmul(nx[i],X[j])))
            row.append(val)
        Q.append(tuple(row))
    return tuple(Q)

Bs=[
 ('diag12',((one,z,z),(z,-one,z),(z,z,z))),
 ('diag21',((arb(2),z,z),(z,-one,z),(z,z,-one))),
 ('mixed',((arb('0.4'),arb('0.7'),arb('-0.2')),(arb('0.7'),arb('-0.1'),arb('0.3')),(arb('-0.2'),arb('0.3'),arb('-0.3')))),
]
rows=[]
for name,B in Bs:
    if not tr(B).contains(0): raise AssertionError(('B trace',name,tr(B)))
    # certify symmetry
    for i in range(3):
        for j in range(3):
            if not (B[i][j]-B[j][i]).contains(0): raise AssertionError(('B symmetry',name,i,j))
    omega=matvec(B,X)
    u=C.vscale(-arb(1)/3,C.cross(X,omega))
    curlerr=C.vadd(C.curl(u),C.vscale(-1,omega)); divu=C.div(u); normal=C.vdot(X,u)
    if not vnorm2(curlerr).contains(0): raise AssertionError(('curl lift',name,vnorm2(curlerr)))
    if not C.savg(C.pmul(divu,divu)).contains(0): raise AssertionError(('div lift',name))
    if not C.savg(C.pmul(normal,normal)).contains(0): raise AssertionError(('tangent lift',name))
    Q0=transaction_Q(omega)
    if not mnorm2(Q0).contains(0): raise AssertionError(('linear transaction null',name,Q0))
    lap=lapv(omega)
    if not vnorm2(lap).contains(0): raise AssertionError(('linear viscous zero mode',name,vnorm2(lap)))
    G=C.bracket(omega,u)
    QG=transaction_Q(G)
    B2tf=mtf(mprod(B,B)); pred=mscale(arb(2)/5,B2tf)
    err=mnorm2(msub(QG,pred))
    if not err.contains(0): raise AssertionError(('self productive Q identity',name,QG,pred,err))
    # For Q_G(r)=r^2 QG on radius r, Hodge screen gives Sdot_v=(5/14)r^2 QG=(1/7)r^2(B^2)_TF.
    Sdot_unit=mscale(arb(5)/14,QG); Spred=mscale(arb(1)/7,B2tf)
    serr=mnorm2(msub(Sdot_unit,Spred))
    if not serr.contains(0): raise AssertionError(('strain source identity',name,serr))
    rows.append({
        'B_case':name,'B_Frobenius_squared':str(mnorm2(B)),
        'linear_vorticity_transaction_Q_squared':str(mnorm2(Q0)),
        'linear_vorticity_laplacian_mean_square':str(vnorm2(lap)),
        'self_Euler_source_mean_square_on_unit_sphere':str(vnorm2(G)),
        'self_Euler_productive_Q':[[str(x) for x in row] for row in QG],
        'two_fifths_B2_TF':[[str(x) for x in row] for row in pred],
        'Q_identity_error_squared':str(err),
        'generated_vortical_strain_rate_coefficient_at_r1':[[str(x) for x in row] for row in Sdot_unit],
        'one_seventh_B2_TF':[[str(x) for x in row] for row in Spred],
    })

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
   'A symmetric-traceless linear vorticity gradient omega=B x is a genuine smooth transaction-null and viscosity-null local mode: its Hodge transaction tensor Q vanishes by parity and Delta omega=0.  Its exact tangent div-curl velocity in the ball is u=-(1/3)x cross Bx. '
   'Nevertheless its self Euler vorticity source G=(omega.grad)u-(u.grad)omega is quadratic and has the exact productive transaction Q_G(r)=(2/5) r^2 (B^2)_TF.  Hodge screening of this regular r^2 source gives generated vortical strain-rate Sdot_v(r)=(1/7) r^2 (B^2)_TF. '
   'Thus transaction-null ancestry is not dynamically inert: the cheapest degree-one poloidal null germ can be a reusable nonlinear catalyst that emits the leading productive degree-two germ while its own leading coefficient is untouched by viscosity. '
   'This kills any argument that treats the transaction-null sector as a purely wasteful positive remainder.  The next question is physical scaling: can such a null catalyst generate the epsilon^-2 maintenance rate of a tiny core while staying finite-energy and low-circulation-Re at its own halo scale?'),
 'rows':rows
},indent=2,allow_nan=False))
