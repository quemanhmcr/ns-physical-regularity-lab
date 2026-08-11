import json, os, math
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import degree6_hodge_servo_core as C
z=C.z; one=C.o; rt2=arb(2).sqrt()
X=({(1,0,0):one},{(0,1,0):one},{(0,0,1):one})
r2=C.padd(C.padd(C.pmul(X[0],X[0]),C.pmul(X[1],X[1])),C.pmul(X[2],X[2]))

def ppow(P,n):
    q={(0,0,0):one}
    for _ in range(n): q=C.pmul(q,P)
    return q

def lin(v):
    q={}
    for j in range(3): q=C.padd(q,C.pscale(v[j],X[j]))
    return q

def Hplane(l,U,V):
    q={}
    for k in range(0,l+1,2):
        q=C.padd(q,C.pscale(((-1)**(k//2))*math.comb(l,k),C.pmul(ppow(U,l-k),ppow(V,k))))
    return q

def matvec(S,V):
    out=[]
    for i in range(3):
        q={}
        for j in range(3): q=C.padd(q,C.pscale(S[i][j],V[j]))
        out.append(q)
    return tuple(out)

def affine_euler(S,V):
    SV=matvec(S,V); Sx=matvec(S,X); adv=[]
    for i in range(3):
        q={}
        for j in range(3): q=C.padd(q,C.pmul(Sx[j],C.pder(V[i],j)))
        adv.append(q)
    return C.vadd(SV,C.vscale(-1,tuple(adv)))

def proj_coeff(T,G):
    den=C.savg(C.vdot(T,T)); num=C.savg(C.vdot(T,G))
    if den.contains(0): raise AssertionError('zero toroidal norm')
    return num/den

def Cscaled(l,coeff,e):
    ep=e**(2*l+1); q=z
    for m,c in enumerate(coeff):
        q += c*(one/arb(2*m+2)-ep/arb(2*m+2*l+3))
    return -arb(l+1)/arb(2*l+1)*q

# Validated stationary productive strain S_*/sigma and eigenframe.
a=arb(3)/(2*rt2)
S=((one,z,a),(z,one,a),(a,a,-2*one))
cplus=(3*rt2-one)/2; cminus=(one+3*rt2)/2
e_mid=(one/rt2,-one/rt2,z)
t=2-rt2; np=(4*(2-rt2)).sqrt(); e_plus=(one/np,one/np,t/np)
u=2+rt2; nm=(4*(2+rt2)).sqrt(); e_minus=(one/nm,one/nm,-u/nm)
frame=[e_plus,e_mid,e_minus]; coords=[lin(v) for v in frame]
pairs=[('plus_mid',0,1),('mid_minus',1,2),('minus_plus',2,0)]

rows=[]
for l in (2,4,6,8):
    lobe_radial=arb(8)/(arb(l+2)*arb(l+4)*arb(l+6))
    for label,i,j in pairs:
        H=Hplane(l,coords[i],coords[j]); grad=tuple(C.pder(H,k) for k in range(3)); T=C.cross(X,grad)
        pq=[]
        for qdeg in (0,2):
            V=tuple(C.pmul(ppow(r2,qdeg//2),v) for v in T)
            pq.append(proj_coeff(T,affine_euler(S,V)))
        alpha=pq[0]; beta=(pq[1]-pq[0])/2
        h=[alpha,-2*alpha-4*beta,alpha+4*beta]; f=[one,-2*one,one]
        for es in ('0.1','1e-2','1e-4','1e-8','1e-12','1e-20'):
            e=arb(es); Craw=Cscaled(l,f,e)
            if Craw.contains(0): raise AssertionError(('screen normalization',l,label,es))
            Amp=one/Craw
            Junit=Amp*Cscaled(l,h,e)
            if Junit.contains(0):
                rows.append({'l':l,'orientation':label,'epsilon':es,'affine_source_blind':True})
                continue
            M1=arb(l+1)*Amp/(e*e)
            sigma_req=-M1/Junit # nu=1
            actor_Re=abs(sigma_req)*e*e
            Gamma_lobe=2*abs(Amp)*(e**l)*lobe_radial
            gap=actor_Re/Gamma_lobe
            scaled_gap=gap*(e**l)
            if not (actor_Re.lower()>0 and Gamma_lobe.lower()>0 and gap.lower()>0):
                raise AssertionError(('positive transaction gap',l,label,es,actor_Re,Gamma_lobe,gap))
            rows.append({
                'l':l,'orientation':label,'epsilon':es,'affine_source_blind':False,
                'maintenance_actor_transaction_Re_sigma_epsilon2_over_nu':str(actor_Re),
                'tiny_lobe_circulation_Re_Gamma_lobe_over_nu':str(Gamma_lobe),
                'actor_transaction_Re_over_lobe_circulation_Re':str(gap),
                'gap_times_epsilon_power_l':str(scaled_gap),
                'sigma_required_times_epsilon_squared_nu1':str(sigma_req*e*e),
                'lobe_Gamma_over_epsilon_power_l':str(Gamma_lobe/(e**l)),
            })

# Each tested orientation is nonblind by module160; enforce same here.
for l in (2,4,6,8):
    for label,_,_ in pairs:
        rr=[r for r in rows if r['l']==l and r['orientation']==label]
        if not rr or any(r['affine_source_blind'] for r in rr):
            raise AssertionError(('unexpected affine blind orientation',l,label))

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
   'Combine the actual affine Euler maintenance gate of module160 with the actual Stokes lobe circulation of module162 for the same compact unit-screened core. '
   'For every tested stationary-strain eigenplane orientation, the maintenance-capable affine actor has transaction Reynolds R_actor=|sigma_required| epsilon^2/nu tending to a finite nonzero angular constant, while the core own material-circulation observable R_lobe=Gamma_lobe/nu scales as epsilon^l and vanishes. '
   'Hence R_actor/R_lobe grows like epsilon^-l; the reported product (R_actor/R_lobe) epsilon^l tends to a finite angular constant. '
   'This is a physical separation between the transaction strength required to maintain the screened core and the circulation carried by one of the core own vorticity lobes. '
   'It is deliberately not identified with a lower bound on independent material circulation stock: through-going winding can amplify transaction circulation relative to one lineage circulation.  The next attack must therefore test whether winding can supply the order-one maintenance transaction Reynolds with vanishing lineage circulation without losing viscous persistence.'),
 'rows':rows
},indent=2,allow_nan=False))
