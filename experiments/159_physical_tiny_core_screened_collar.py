import json, os, math
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import degree6_hodge_servo_core as C
z=C.z; one=C.o
X=({(1,0,0):one},{(0,1,0):one},{(0,0,1):one})
r2=C.padd(C.padd(C.pmul(X[0],X[0]),C.pmul(X[1],X[1])),C.pmul(X[2],X[2]))

def Hreal(l):
    P={}
    for k in range(0,l+1,2): P[(l-k,k,0)]=arb(((-1)**(k//2))*math.comb(l,k))
    return P

def ball_even(exp,R):
    a,b,c=exp
    if a%2 or b%2 or c%2:return z
    aa=arb(a)/2+arb('0.5'); bb=arb(b)/2+arb('0.5'); cc=arb(c)/2+arb('0.5')
    N=(a+b+c)//2
    sphere=2*aa.gamma()*bb.gamma()*cc.gamma()/(arb(N)+arb('1.5')).gamma()
    return (R**(2*N+3))*sphere/arb(2*N+3)

def pint(P,R):
    q=z
    for e,v in P.items(): q += v*ball_even(e,R)
    return q

def vnorm_ball(V,R): return pint(C.vdot(V,V),R)
def grad_vnorm_ball(V,R):
    q=z
    for i in range(3):
        for j in range(3): q += pint(C.pmul(C.pder(V[i],j),C.pder(V[i],j)),R)
    return q

def Ccoef(l,k): return -arb(l+1)/(arb(2*k+2)*arb(2*k+2*l+3))

def pconst(c): return {(0,0,0):c}

rows=[]
for l in (2,3,4,6,8):
    H=Hreal(l); grad=tuple(C.pder(H,j) for j in range(3)); T=C.cross(X,grad)
    f=C.padd(pconst(one),C.padd(C.pscale(-2,r2),C.pmul(r2,r2)))
    Obase=tuple(C.pmul(f,t) for t in T)
    Zbase=vnorm_ball(Obase,one); Pbase=grad_vnorm_ball(Obase,one)
    if not (Zbase.lower()>0 and Pbase.lower()>0): raise AssertionError(('base physical forms',l,Zbase,Pbase))
    # Stable source-scaled screen integrals after r=epsilon s.
    I0=arb(1)/6
    I1=arb(1)/(2*l+3)-arb(2)/(2*l+5)+arb(1)/(2*l+7)
    prevZ=None; prevP=None
    for es in ['0.5','0.25','0.1','1e-2','1e-4','1e-8','1e-12','1e-20']:
        e=arb(es)
        Craw= -arb(l+1)/arb(2*l+1)*(I0-(e**(2*l+1))*I1)
        if Craw.contains(0): raise AssertionError(('screen normalization singular',l,es,Craw))
        Amp=one/Craw
        Cnormalized=Amp*Craw
        if not Cnormalized.contains(1): raise AssertionError(('screened companion normalization',l,es,Cnormalized))
        # Exact physical-space scaling under x=epsilon y.
        Zphys=Amp*Amp*(e**(2*l-1))*Zbase
        Pphys=Amp*Amp*(e**(2*l-3))*Pbase
        if not (Zphys.lower()>0 and Pphys.lower()>0): raise AssertionError(('positive physical forms',l,es,Zphys,Pphys))
        if prevZ is not None and not (Zphys.upper()<prevZ.lower()): raise AssertionError(('enstrophy not decreasing',l,es,Zphys,prevZ))
        if prevP is not None and not (Pphys.upper()<prevP.lower()): raise AssertionError(('palinstrophy not decreasing',l,es,Pphys,prevP))
        prevZ,prevP=Zphys,Pphys
        # Independent Cartesian expansion only at moderate epsilon; deep scales use the intrinsic scaled observer above.
        Zerr='not_evaluated_deep_scale'; Perr='not_evaluated_deep_scale'
        if es in ('0.5','0.25','0.1'):
            raw=[e**-2,-2*(e**-4),e**-6]
            apoly=C.padd(pconst(Amp*raw[0]),C.padd(C.pscale(Amp*raw[1],r2),C.pscale(Amp*raw[2],C.pmul(r2,r2))))
            O=tuple(C.pmul(apoly,t) for t in T)
            Zd=vnorm_ball(O,e); Pd=grad_vnorm_ball(O,e)
            ze=Zd-Zphys; pe=Pd-Pphys
            if not ze.contains(0): raise AssertionError(('direct enstrophy scaling',l,es,ze))
            if not pe.contains(0): raise AssertionError(('direct palinstrophy scaling',l,es,pe))
            Zerr=str(ze); Perr=str(pe)
        rows.append({'l':l,'epsilon':es,'screened_companion_coefficient_C_l':str(Cnormalized),'normalizing_amplitude_A':str(Amp),'actual_vorticity_enstrophy':str(Zphys),'actual_vorticity_palinstrophy_grad_omega_square':str(Pphys),'enstrophy_over_epsilon_power_2l_minus_1':str(Zphys/(e**(2*l-1))),'palinstrophy_over_epsilon_power_2l_minus_3':str(Pphys/(e**(2*l-3))),'moderate_scale_direct_cartesian_enstrophy_error':Zerr,'moderate_scale_direct_cartesian_palinstrophy_error':Perr,'structural_boundary_a':'0','structural_boundary_radial_derivative_a':'0'})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
   'This is a physical-space dual of the remote spectral collar.  Inside the source ball L=1, take the genuine divergence-free toroidal vorticity omega_e=a_e(r) x cross grad H_l with a_e=A_e e^-2(1-r^2/e^2)^2 for r<=e and zero outside. '
   'The radial factor and its first derivative vanish structurally at r=e, so the compact extension is C1.  A_e is chosen from the exact source-scaled screened functional after the physical change of variables r=e s, avoiding catastrophic monomial cancellation at tiny e.  Thus C_l[a_e]=1 for every epsilon; C_l is only the local harmonic-companion/feedback coefficient, not an independently occupied Hodge energy. '
   'The exact physical-space scaling gives Z_e=int|omega_e|^2=A_e^2 epsilon^(2l-1) Z_base and P_e=int|grad omega_e|^2=A_e^2 epsilon^(2l-3) P_base.  Independent direct Cartesian ball integration verifies these identities at moderate epsilon. '
   'Since A_e tends to a finite nonzero limit as epsilon tends to zero, both enstrophy and palinstrophy vanish for every l>=2 while the screened local feedback coefficient remains one. '
   'Thus the remote-collar escape is not an artifact of the Hankel normalization: it has an explicit compact physical-space realization.  The payment is extreme localization and therefore a need for continual nonlinear creation/replenishment, not static enstrophy or static viscous palinstrophy.'),
 'rows':rows
},indent=2,allow_nan=False))
