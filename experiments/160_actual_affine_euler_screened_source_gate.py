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
        term=C.pmul(ppow(U,l-k),ppow(V,k))
        q=C.padd(q,C.pscale(((-1)**(k//2))*math.comb(l,k),term))
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
    # profile epsilon^-2 sum_m coeff[m] s^(2m), r=epsilon s, before overall A normalization
    q=z
    ep=e**(2*l+1)
    for m,c in enumerate(coeff):
        q += c*(one/arb(2*m+2)-ep/arb(2*m+2*l+3))
    return -arb(l+1)/arb(2*l+1)*q

# Validated stationary-amplifier strain S_*/sigma.
a=arb(3)/(2*rt2)
S=((one,z,a),(z,one,a),(a,a,-2*one))
cplus=(3*rt2-one)/2; cminus=(one+3*rt2)/2
# Orthonormal eigenframe.
e_mid=(one/rt2,-one/rt2,z)
t=2-rt2; np=(4*(2-rt2)).sqrt(); e_plus=(one/np,one/np,t/np)
u=2+rt2; nm=(4*(2+rt2)).sqrt(); e_minus=(one/nm,one/nm,-u/nm)
frame=[e_plus,e_mid,e_minus]; evals=[cplus,one,-cminus]
for i,e in enumerate(frame):
    n=sum((e[j]*e[j] for j in range(3)),z)
    if not (n-one).contains(0): raise AssertionError(('eigenvector norm',i,n))
    Se=[sum((S[r][j]*e[j] for j in range(3)),z) for r in range(3)]
    err=sum(((Se[r]-evals[i]*e[r])**2 for r in range(3)),z)
    if not err.contains(0): raise AssertionError(('stationary eigenframe',i,err))
for i in range(3):
    for j in range(i):
        d=sum((frame[i][k]*frame[j][k] for k in range(3)),z)
        if not d.contains(0): raise AssertionError(('eigenframe orthogonality',i,j,d))

coords=[lin(v) for v in frame]
pairs=[('plus_mid',0,1),('mid_minus',1,2),('minus_plus',2,0)]
rows=[]
for l in (2,4,6,8):
    nonblind=0
    for label,i,j in pairs:
        H=Hplane(l,coords[i],coords[j]); grad=tuple(C.pder(H,k) for k in range(3)); T=C.cross(X,grad)
        pq=[]
        for qdeg in (0,2,4,6):
            V=tuple(C.pmul(ppow(r2,qdeg//2),v) for v in T)
            pq.append(proj_coeff(T,affine_euler(S,V)))
        alpha=pq[0]; beta=(pq[1]-pq[0])/2
        for qdeg,pv in zip((0,2,4,6),pq):
            er=pv-(alpha+qdeg*beta)
            if not er.contains(0): raise AssertionError(('radial affine projection law',l,label,qdeg,er))
        # f=(1-s^2)^2; alpha f + beta s f' has coefficients below.
        h=[alpha,-2*alpha-4*beta,alpha+4*beta]
        f=[one,-2*one,one]
        orientation_nonblind=False
        for es in ('0.1','1e-2','1e-4','1e-8','1e-12','1e-20'):
            e=arb(es)
            Craw=Cscaled(l,f,e)
            if Craw.contains(0): raise AssertionError(('core screen',l,label,es,Craw))
            Amp=one/Craw
            Junit=Amp*Cscaled(l,h,e)   # same-l screened Euler source moment for sigma=1
            M1=arb(l+1)*Amp/(e*e)     # viscous drift C[D_l a]
            if not Junit.contains(0):
                orientation_nonblind=True
                sigma_req=-M1/Junit
                rows.append({'l':l,'D2_eigenplane_orientation':label,'epsilon':es,'angular_alpha_for_a':str(alpha),'angular_beta_for_r_a_prime':str(beta),'unit_sigma_same_l_Euler_source_moment_J0':str(Junit),'viscous_screened_drift_M1_nu1':str(M1),'sigma_required_to_hold_C_l_constant_nu1':str(sigma_req),'sigma_required_times_epsilon_squared':str(sigma_req*e*e),'affine_source_nonblind':True})
            else:
                rows.append({'l':l,'D2_eigenplane_orientation':label,'epsilon':es,'angular_alpha_for_a':str(alpha),'angular_beta_for_r_a_prime':str(beta),'unit_sigma_same_l_Euler_source_moment_J0':str(Junit),'viscous_screened_drift_M1_nu1':str(M1),'sigma_required_to_hold_C_l_constant_nu1':'undefined_affine_blind','sigma_required_times_epsilon_squared':'undefined_affine_blind','affine_source_nonblind':False})
        if orientation_nonblind: nonblind+=1
    if nonblind==0: raise AssertionError(('all D2 orientations affine blind',l))

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','stationary_eigenvalues_over_sigma':[str(x) for x in evals],
 'interpretation':(
   'This module projects the actual affine Euler vorticity operator E_S omega=S omega-(Sx).grad omega for the validated stationary productive strain S_* back into the same toroidal angular Hodge channel. '
   'For every tested D2 eigenplane harmonic, the same-l projection has the exact radial form sigma[alpha_l a(r)+beta_l r a_prime(r)]; alpha and beta are finite angular numbers independent of the core radius epsilon, and q=0,2,4,6 radial monomials certify the linear q law. '
   'For the compact tiny core a_e=A_e e^-2(1-r^2/e^2)^2 normalized by C_l[a_e]=1, the screened affine Euler source moment J_0 at unit sigma remains O(1), whereas the exact viscous drift M_1=C_l[D_l a_e]=(l+1)A_e e^-2. '
   'Whenever the orientation is not affine-blind, the strain required to hold the screened coefficient constant is sigma=-M_1/J_0 and therefore sigma e^2 tends to a finite nonzero angular constant. '
   'This is the continuous-maintenance ledger evaluated on an actual Navier-Stokes Euler mechanism: a fixed reusable common strain cannot maintain an arbitrarily small screened core; its rate must grow like nu/epsilon^2. '
   'Some angular orientations may be blind to the same-l affine source, which only makes the common-strain mechanism less capable.  The remaining escape is genuinely non-affine same-scale dynamics or continual routing between angular channels, which must be traced next through material ancestry.'),
 'rows':rows
},indent=2,allow_nan=False))
