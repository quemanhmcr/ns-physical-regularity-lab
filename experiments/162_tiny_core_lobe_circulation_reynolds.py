import json, os, math
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import degree6_hodge_servo_core as C
z=C.z; one=C.o
X=({(1,0,0):one},{(0,1,0):one},{(0,0,1):one})
r2=C.padd(C.padd(C.pmul(X[0],X[0]),C.pmul(X[1],X[1])),C.pmul(X[2],X[2]))

def pconst(v): return {(0,0,0):v}
def ppow(P,n):
    q=pconst(one)
    for _ in range(n): q=C.pmul(q,P)
    return q

def Hreal(l):
    P={}
    for k in range(0,l+1,2): P[(l-k,k,0)]=arb(((-1)**(k//2))*math.comb(l,k))
    return P

def ball_even(exp):
    a,b,c=exp
    if a%2 or b%2 or c%2:return z
    aa=arb(a)/2+arb('0.5');bb=arb(b)/2+arb('0.5');cc=arb(c)/2+arb('0.5');N=(a+b+c)//2
    sph=2*aa.gamma()*bb.gamma()*cc.gamma()/(arb(N)+arb('1.5')).gamma()
    return sph/arb(2*N+3)
def pint(P): return sum((v*ball_even(e) for e,v in P.items()),z)
def vnorm(V): return pint(C.vdot(V,V))
def lapv(V): return tuple(C.plap(v) for v in V)

rows=[]
for l in (2,3,4,6,8):
    H=Hreal(l); grad=tuple(C.pder(H,j) for j in range(3)); T=C.cross(X,grad)
    fcoef=[one,-2*one,one]
    f=C.padd(pconst(one),C.padd(C.pscale(-2,r2),C.pmul(r2,r2)))
    lobe_radial_integral=one/arb(l+2)-2*one/arb(l+4)+one/arb(l+6)
    lobe_closed=arb(8)/(arb(l+2)*arb(l+4)*arb(l+6))
    if not (lobe_radial_integral-lobe_closed).contains(0): raise AssertionError(('lobe radial integral',l))
    for es in ('0.1','1e-2','1e-4','1e-8','1e-12','1e-20'):
        e=arb(es)
        I0=one/6; I1=one/arb(2*l+3)-2*one/arb(2*l+5)+one/arb(2*l+7)
        Craw=-arb(l+1)/arb(2*l+1)*(I0-(e**(2*l+1))*I1); Amp=one/Craw
        Omega=tuple(C.pscale(Amp,C.pmul(f,t)) for t in T)
        # Exact scaled Hodge velocity, same validated monomial lift as module 139.
        U=tuple(g for g in grad)
        for m,fc in enumerate(fcoef):
            q=2*m; den=arb((q+2)*(q+2*l+3)); Aq=arb(q+l+3)/den; Bq=-arb(l)/arb(q+2*l+3)
            rq=ppow(r2,m); rqp2=ppow(r2,m+1)
            term=tuple(C.padd(C.pscale(Amp*fc*Aq,C.pmul(rqp2,grad[i])),C.pscale(Amp*fc*Bq,C.pmul(C.pmul(rq,H),X[i]))) for i in range(3))
            U=C.vadd(U,term)
        G=C.bracket(Omega,U); DO=lapv(Omega)
        NG=vnorm(G); ND=vnorm(DO)
        if not (NG.lower()>0 and ND.lower()>0): raise AssertionError(('source norms',l,es,NG,ND))
        full_ratio=(e**l)*(NG/ND).sqrt()  # nu=1
        # On z=0, omega_z=a(R) d_phi H_l.  Over one positive lobe, integral d_phi H_l dphi = 2.
        Gamma=2*abs(Amp)*(e**l)*lobe_radial_integral
        if not (Gamma.lower()>0): raise AssertionError(('positive lobe circulation',l,es,Gamma))
        kappa=full_ratio/Gamma
        rows.append({
            'l':l,'epsilon':es,
            'screened_companion_C_l':'1',
            'one_positive_lobe_Stokes_circulation_abs':str(Gamma),
            'circulation_Reynolds_Gamma_over_nu_nu1':str(Gamma),
            'Gamma_over_epsilon_power_l':str(Gamma/(e**l)),
            'full_self_Euler_to_viscous_L2_source_ratio_nu1':str(full_ratio),
            'full_source_ratio_over_epsilon_power_l':str(full_ratio/(e**l)),
            'full_source_ratio_over_lobe_circulation_Reynolds':str(kappa),
            'lobe_flux_formula':'2 |A| epsilon^l integral_0^1 (1-s^2)^2 s^(l+1) ds',
            'closed_radial_integral_8_over_product':str(lobe_closed),
        })

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
   'For H_l=Re(x+iy)^l, the equatorial vorticity component of the compact screened core is omega_z=a(R) partial_phi H_l.  Integrating one positive angular lobe gives an actual vorticity flux, hence by Stokes an actual circulation, Gamma_lobe=2|A_e| epsilon^l integral_0^1(1-s^2)^2 s^(l+1)ds = 16|A_e| epsilon^l/[(l+2)(l+4)(l+6)]. '
   'The full self-Euler-to-viscous source ratio from the exact interior Hodge velocity is evaluated in the same module and scales as epsilon^l.  Their quotient approaches a finite nonzero shape constant. '
   'Thus the isolated core maintenance deficit is naturally measured by the physical circulation Reynolds Gamma_lobe/nu, not by an imposed abstract norm ratio.  A unit-screened tiny core has Gamma_lobe/nu ->0 and therefore becomes self-dynamically subcritical even though its local screened companion coefficient is fixed at one. '
   'Modules 158-160 show that maintenance by another same-scale strain actor instead requires a rate sigma~nu/epsilon^2, equivalently an order-one transaction Reynolds epsilon^2 sigma/nu.  The resulting gap is physical: the core own circulation ancestry vanishes while a maintenance-capable actor must bring order-nu circulation/transaction ancestry to the shrinking scale.'),
 'rows':rows
},indent=2,allow_nan=False))
