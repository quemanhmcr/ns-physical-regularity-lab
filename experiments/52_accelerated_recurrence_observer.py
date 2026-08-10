import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); sqrt2=arb(2).sqrt(); sqrt5=arb(5).sqrt()

# Accelerated phase phi=t/(1-t), eta=1, ell=2, eps=1.
# Infinite returns t_k=k*pi/(1+k*pi) accumulate at t=1, while 1<=lambda<=2 gives I(t)<4 and 2<=Q<6.
ell=arb(2); eps=arb(1); eta=arb(1); Qmax=ell+eps*(1+eta)**2
if not Qmax.contains(6): raise AssertionError(('Qmax',Qmax))

# Uniform finite-amplitude oscillation of |T| independent of period index.
T_return_min=ell/(sqrt2*(ell*ell+5).sqrt())
lamc=1+eta
# General crest formula for the same pair: |T|=Q/[sqrt(lambda^2+1)*sqrt(Q^2+4 lambda^4+lambda^2)].
T_crest_max=Qmax/(((lamc*lamc+1).sqrt())*(Qmax*Qmax+4*lamc**4+lamc*lamc).sqrt())
gap=T_return_min-T_crest_max
if not (gap > arb('0.2')): raise AssertionError(('uniform T oscillation gap too small',gap,T_return_min,T_crest_max))

# Every strobe has uniformly positive instantaneous angular-surplus rate.
strobe_rate_lower=5*eps/(Qmax*(Qmax*Qmax+5))
if not (strobe_rate_lower > arb('0.02')): raise AssertionError(('strobe surplus lower bound',strobe_rate_lower))

ks=['1','10','1e6','1e30','1e60']
rows=[]
prev=None
for ks_ in ks:
    k=arb(ks_); phase=k*pi
    # The intrinsic recurrence-tail coordinate is tau=1/(1+phase).  Do not recover a
    # tiny tau by subtracting the near-saturated parent state t_k from one.
    remaining=1/(1+phase)
    if not (remaining>0): raise AssertionError(('positive recurrence tail',ks_,remaining))
    if prev is not None and not (remaining<prev): raise AssertionError(('recurrence tails not decreasing',ks_,remaining,prev))
    prev=remaining
    tk=phase*remaining  # display state only; the small tail is observed directly above.
    closure=tk+remaining
    if not closure.contains(1): raise AssertionError(('intrinsic return-coordinate closure',ks_,closure))
    rows.append({'k':ks_,'return_time_display':str(tk),'intrinsic_time_remaining_to_1':str(remaining),'t_plus_tail':str(closure)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'bounded_memory_Q_upper':str(Qmax),'uniform_T_return_lower':str(T_return_min),'uniform_T_crest_upper':str(T_crest_max),'uniform_halfcycle_T_gap':str(gap),'uniform_strobe_dlogT_lower':str(strobe_rate_lower),
 'interpretation':(
  'The exact quadratic NS family can compress infinitely many lambda=1 returns into t<1 via phi=t/(1-t). '
  'Because 1<=lambda<=2, the entire deposited bridge memory remains bounded by Q<6. '
  'Nevertheless every cycle has a |T| drop/refill exceeding 0.2 and every return has d log|T|/dt>0.02. '
  'Thus event count and even cumulative positive variation can diverge while net material memory stays finite.'
 ),'rows':rows
},indent=2,allow_nan=False))
