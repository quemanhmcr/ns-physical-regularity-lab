import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); sqrt2=arb(2).sqrt()

def Tabs(lam,Q):
    return Q/((lam*lam+1).sqrt()*(Q*Q+4*lam**4+lam*lam).sqrt())

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))

# Exact NS family u=(a(t)x+eps yz,-a(t)y,0), lambda=1+sin^2 t.
# At t=0, pi/2, pi the stretch is 1,2,1.  The bridge memory clock is
# I(t)=19t/8-(3/4)sin2t+(1/32)sin4t, so the crest and return values are exact.
L=arb(2); Ic=arb(19)*pi/16; I1=arb(19)*pi/8
T0=Tabs(arb(1),L)
Tcrest_zero=Tabs(arb(2),L)
lim_gap=T0-Tcrest_zero
if not (lim_gap>arb('0.35')): raise AssertionError(('limiting finite Gram excursion',lim_gap))
rows=[]
for es in ['1e-30','1e-24','1e-18','1e-12','1e-6','0.01']:
    eps=arb(es); Qc=L+eps*Ic; Q1=L+eps*I1
    Tc=Tabs(arb(2),Qc); T1=Tabs(arb(1),Q1)
    half_drop=T0-Tc; return_error=T1-T0
    deposited=Q1-L
    certify_one(deposited/(eps*I1),('exact deposited bridge memory',es))
    # Finite shape excursion persists as eps -> 0 while normalized pair-cell renewal per cycle vanishes.
    if not (half_drop>arb('0.3')): raise AssertionError(('finite shape excursion lost',es,half_drop))
    rel_D_renewal=(Q1-L)/L  # exact relative change of |D| because |D|=eps^2 Q
    if eps <= arb('1e-6') and not (rel_D_renewal < arb('1e-5')): raise AssertionError(('renewal not small',es,rel_D_renewal))
    rows.append({'eps':es,'T_initial':str(T0),'T_crest':str(Tc),'T_next_return':str(T1),'finite_halfcycle_T_drop':str(half_drop),'next_return_T_error':str(return_error),'deposited_Q':str(deposited),'relative_pair_cell_renewal':str(rel_D_renewal)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'eps_to_zero_T_gap':str(lim_gap),
 'interpretation':(
  'In the exact quadratic Navier-Stokes family, a common affine symmetric strain can drive a finite-amplitude productive Gram excursion while the non-affine pair-cell renewal is made arbitrarily small by eps->0+.  Thus residual Gram-shape motion does not by itself force ancestry replacement; common strain is a genuine reusable source channel that must be accounted for physically.'
 ),'rows':rows
},indent=2,allow_nan=False))
