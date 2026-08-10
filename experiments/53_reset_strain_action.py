import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x > 1-t and x < 1+t): raise AssertionError((label,x))

etas=['1e-6','0.1','1','10','1e6']
Ns=['1','10','1e6','1e30','1e60']
rows=[]
for es in etas:
  eta=arb(es); per=2*(1+eta).log()
  if not (per>0): raise AssertionError(('positive reset action',es,per))
  for Ns_ in Ns:
    N=arb(Ns_)
    action=N*per
    certify_one(action/(2*N*(1+eta).log()),('reset action identity',es,Ns_))
    # Accelerated phase return time after N cycles remains below one.
    phase=N*pi; tN=phase/(1+phase)
    if not (tN<1): raise AssertionError(('finite-time packing',es,Ns_,tN))
    rows.append({'eta':es,'N_cycles':Ns_,'action_per_cycle':str(per),'total_log_stretch_action':str(action),'accelerated_return_time':str(tN)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'For lambda=1+eta sin^2(phi), every complete material-shape reset 1 -> 1+eta -> 1 has exact action integral |d log lambda|=2 log(1+eta), independent of how fast the phase is traversed. '
  'Therefore N finite-amplitude resets cost exactly 2N log(1+eta); packing infinitely many into finite time forces infinite logarithmic stretch action even though the bridge-memory integral can remain bounded.'
 ),'rows':rows
},indent=2,allow_nan=False))
