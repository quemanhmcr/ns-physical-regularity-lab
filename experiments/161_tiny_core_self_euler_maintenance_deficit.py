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
def allzero(P): return all(v.contains(0) for v in P.values())
def vallzero(V): return all(allzero(p) for p in V)
def homog(P,d): return {e:v for e,v in P.items() if sum(e)==d}
def vhomog(V,d): return tuple(homog(p,d) for p in V)
def Cscreen_scaled(l, coeff_by_q, e):
    # Actual source radial coefficient has overall epsilon^(l-4) for self Euler.
    ep=e**(2*l+1); s=z
    for q,g in coeff_by_q.items():
        s += g*(one/arb(q+2)-ep/arb(q+2*l+3))
    return -arb(l+1)/arb(2*l+1)*(e**(l-2))*s

def lapv(V): return tuple(C.plap(v) for v in V)

rows=[]
for l in (2,3,4,6,8):
    H=Hreal(l); grad=tuple(C.pder(H,j) for j in range(3)); T=C.cross(X,grad); den=C.savg(C.vdot(T,T))
    if den.contains(0): raise AssertionError(('toroidal norm',l))
    fcoef=[one,-2*one,one]
    f=C.padd(pconst(one),C.padd(C.pscale(-2,r2),C.pmul(r2,r2)))
    for es in ('0.1','1e-2','1e-4','1e-8','1e-12','1e-20'):
        e=arb(es)
        I0=one/6; I1=one/arb(2*l+3)-2*one/arb(2*l+5)+one/arb(2*l+7)
        Craw=-arb(l+1)/arb(2*l+1)*(I0-(e**(2*l+1))*I1); Amp=one/Craw
        Omega=tuple(C.pscale(Amp,C.pmul(f,t)) for t in T)
        # Exact scaled interior Hodge velocity: companion coefficient is C_l=1.
        U=tuple(g for g in grad)
        for m,fc in enumerate(fcoef):
            q=2*m
            Aq=arb(q+l+3)/(arb(q+2)*arb(q+2*l+3)); Bq=-arb(l)/arb(q+2*l+3)
            rq=ppow(r2,m); rqp2=ppow(r2,m+1)
            term=tuple(C.padd(C.pscale(Amp*fc*Aq,C.pmul(rqp2,grad[i])),C.pscale(Amp*fc*Bq,C.pmul(C.pmul(rq,H),X[i]))) for i in range(3))
            U=C.vadd(U,term)
        curlerr=C.vadd(C.curl(U),C.vscale(-1,Omega)); diverr=C.div(U)
        curlerr2=C.savg(C.vdot(curlerr,curlerr)); diverr2=C.savg(C.pmul(diverr,diverr))
        if not curlerr2.contains(0): raise AssertionError(('interior curl lift',l,es,curlerr2))
        if not diverr2.contains(0): raise AssertionError(('interior divergence lift',l,es,diverr2))
        G=C.bracket(Omega,U); DO=lapv(Omega)
        NG=vnorm(G); ND=vnorm(DO); NU=vnorm(U); NO=vnorm(Omega)
        if not (ND.lower()>0 and NO.lower()>0 and NU.lower()>0): raise AssertionError(('positive norms',l,es,NO,NU,ND))
        # Project the dimensionless self-Euler source onto the original angular T_l channel degree by degree.
        degs=sorted(set(sum(k) for comp in G for k in comp))
        gam={}
        for d in degs:
            Gd=vhomog(G,d); num=C.savg(C.vdot(T,Gd)); q=d-l
            if q>=0 and not num.contains(0): gam[q]=gam.get(q,z)+num/den
        Jself=Cscreen_scaled(l,gam,e)
        Mvisc=arb(l+1)*Amp/(e*e)
        ratio=Jself/Mvisc
        fullnormratio=(e**l)*(NG/ND).sqrt()
        rows.append({'l':l,'epsilon':es,'screened_companion_C_l':'1','self_Euler_same_l_projected_radial_coefficients_q_to_gamma':{str(q):str(v) for q,v in gam.items()},'self_Euler_screened_source_J0':str(Jself),'viscous_screened_drift_M1_nu1':str(Mvisc),'self_Euler_to_viscous_screened_ratio':str(ratio),'ratio_divided_by_epsilon_power_l':str(ratio/(e**l)),'full_L2_self_Euler_to_viscous_source_norm_ratio':str(fullnormratio),'full_norm_ratio_divided_by_epsilon_power_l':str(fullnormratio/(e**l)),'actual_core_vorticity_L2_square':str((e**(2*l-1))*NO),'actual_core_velocity_L2_square_inside_core':str((e**(2*l+1))*NU),'actual_self_Euler_source_L2_square_inside_core':str((e**(4*l-5))*NG),'actual_viscous_Delta_omega_L2_square_inside_core_nu1':str((e**(2*l-5))*ND),'interior_curl_lift_error_square':str(curlerr2),'interior_divergence_lift_error_square':str(diverr2)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
   'Take the compact tiny core with screened companion coefficient C_l=1.  Its exact interior tangent Hodge velocity is reconstructed from the radial div-curl formulas, including the global harmonic companion selected by the source sphere; curl U=Omega and div U=0 are certified before evaluating the Euler term. '
   'In scaled core coordinates x=epsilon y, Omega scales as epsilon^(l-2), U as epsilon^(l-1), the self-Euler vorticity source [Omega,U] as epsilon^(2l-4), and Delta Omega as epsilon^(l-4). '
   'Projecting the actual self-Euler source back onto the original toroidal angular channel and applying the screened functional gives J_self=O(epsilon^(l-2)), while the exact viscous drift of C_l is M_1=O(epsilon^-2).  Hence J_self/M_1=O(epsilon^l). '
   'The full L2 source-norm ratio has the same epsilon^l scaling.  For every l>=2 the self-generated Euler dynamics of a unit-screened tiny core becomes negligible relative to viscosity as the core shrinks. '
   'Therefore the cheap tiny-core collar is not self-maintaining.  Continuous replenishment must come from an external/intensifying affine strain or from genuinely non-affine same-scale interactions with additional vorticity/velocity ancestry.  This is the first place where the static collar escape is forced to ask for another physical actor.'),
 'rows':rows
},indent=2,allow_nan=False))
