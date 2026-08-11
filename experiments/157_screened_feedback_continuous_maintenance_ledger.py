import json, os
from fractions import Fraction as Q
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def Ccoef(l,k): return -Q(l+1,(2*k+2)*(2*k+2*l+3))
def Dlam(l,k): return Q(0) if k==0 else Q(2*k*(2*k+2*l+1))
def D(l,c):
    if len(c)<=1: return [Q(0)]
    out=[Q(0)]*(len(c)-1)
    for k in range(1,len(c)): out[k-1]+=c[k]*Dlam(l,k)
    return out
def Dq(l,c,q):
    r=c[:]
    for _ in range(q): r=D(l,r)
    return r
def C(l,c): return sum((x*Ccoef(l,k) for k,x in enumerate(c)),Q(0))
def add(a,b):
    n=max(len(a),len(b)); return [(a[i] if i<len(a) else Q(0))+(b[i] if i<len(b) else Q(0)) for i in range(n)]
def scale(a,s): return [s*x for x in a]
def a0(c): return c[0] if c else Q(0)
def a1(c): return sum(c,Q(0))
def Z(l,c): return sum((c[i]*c[j]*Q(1,2*l+2*i+2*j+3) for i in range(len(c)) for j in range(len(c))),Q(0))
def A(q): return arb(q.numerator)/q.denominator

def delayed(l,n):
    # Unique degree-2n polynomial with top coefficient 1 and C[D^q a]=0 for q<n.
    c=[Q(0)]*(n+1); c[n]=Q(1)
    for q in range(n-2,-1,-1):
        # At level q, only coefficients k>=q+1 survive after D^q.
        tail=sum((c[k]*C(l,Dq(l,[Q(0)]*k+[Q(1)],q)) for k in range(q+2,n+1)),Q(0))
        basis=C(l,Dq(l,[Q(0)]*(q+1)+[Q(1)],q))
        c[q+1]=-tail/basis
    tail=sum((c[k]*Ccoef(l,k) for k in range(1,n+1)),Q(0))
    c[0]=-tail/Ccoef(l,0)
    return c

rows=[]
profiles=[
    ([Q(3,7),Q(-5,3),Q(11,13),Q(-2,5),Q(7,19)],[Q(-7,11),Q(4,9),Q(5,17),Q(-3,23)]),
    ([Q(10**20),Q(-3,10**10),Q(7,10**30),Q(-1,10**40),Q(9,10**50)],[Q(-10**15),Q(2,10**7),Q(-5,10**25),Q(1,10**35)]),
]
for l in (2,4,6,8):
    for pi,(a,g) in enumerate(profiles):
        adot=add(g,D(l,a))
        for q in range(5):
            aq=Dq(l,a,q); gq=Dq(l,g,q); dqadot=Dq(l,adot,q)
            Mq=C(l,aq); Jq=C(l,gq); Mn=C(l,D(l,aq))
            chain=C(l,dqadot)-(Jq+Mn)
            endpoint=Mn-Q(l+1)*(a0(aq)-a1(aq))
            if chain!=0: raise AssertionError(('continuous chain',l,pi,q,chain))
            if endpoint!=0: raise AssertionError(('endpoint chain',l,pi,q,endpoint))
            rows.append({'l':l,'profile_case':pi,'q':q,'M_q':str(A(Mq)),'J_q':str(A(Jq)),'M_q_plus_1':str(A(Mn)),'continuous_chain_error':'0','endpoint_identity_error':'0'})

servo=[]
for l in (2,4,6,8):
    for n in range(1,9):
        a=delayed(l,n)
        Ms=[C(l,Dq(l,a,q)) for q in range(n+1)]
        if any(Ms[q]!=0 for q in range(n)) or Ms[n]==0: raise AssertionError(('delayed moments',l,n,Ms))
        g=scale(D(l,a),Q(-1))
        if any(add(g,D(l,a))): raise AssertionError(('freeze source',l,n))
        Js=[C(l,Dq(l,g,q)) for q in range(n)]
        for q in range(n):
            if Js[q]!=-Ms[q+1]: raise AssertionError(('servo',l,n,q,Js[q],Ms[q+1]))
        only_top=all(Js[q]==0 for q in range(n-1)) and Js[-1]!=0
        if not only_top: raise AssertionError(('top-only initial servo',l,n,Js))
        servo.append({'l':l,'hidden_moment_count_n':n,'M_0_through_M_n':[str(A(x)) for x in Ms],'required_J_0_through_J_n_minus_1':[str(A(x)) for x in Js],'only_top_source_moment_nonzero_initially':True,'state_radial_enstrophy_factor':str(A(Z(l,a))),'artificial_freeze_source_radial_L2_factor':str(A(Z(l,g)))})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','generic_chain_cases':len(rows),'delayed_servo_cases':len(servo),
 'interpretation':(
   'In source-scaled time, write one projected toroidal radial channel as partial_tau a_l = g_l + D_l a_l. '
   'Here g_l is the non-viscous source; in actual Navier-Stokes it must be the projection of Euler transport/stretching and angular-channel couplings, not a free forcing. '
   'For M_q=C_l[D_l^q a_l] and J_q=C_l[D_l^q g_l], the exact continuous ledger is dot M_q=J_q+M_(q+1). '
   'The capacitary endpoint law gives M_(q+1)=(l+1)[D_l^q a_l(0)-D_l^q a_l(1)]. '
   'Thus keeping M_0 through M_(n-1) identically zero on a time interval requires the servo J_q=-M_(q+1) at every time. '
   'Maximally delayed polynomial states initially need only the top source moment J_(n-1); the artificial calibration g=-D_l a freezes the whole profile and realizes the chain exactly. '
   'That artificial source is not promoted as an NS mechanism.  It is a microscope showing precisely what the true projected Euler term would have to supply continuously. '
   'Since modules 155-156 make already-created remote collars arbitrarily cheap in static enstrophy and palinstrophy, the remaining physical bottleneck is nonlinear creation/routing and material ancestry of this continuous source ledger.'),
 'generic_rows':rows,'delayed_servo_rows':servo
},indent=2,allow_nan=False))
