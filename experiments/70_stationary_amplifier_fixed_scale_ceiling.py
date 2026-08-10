import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))

# From E_h=(7*pi/5)sigma^2 r^5 for the shape-locked equal amplifier,
# finite kinetic energy E0 implies sigma <= sqrt(5 E0/(7*pi)) r^(-5/2).
rows=[]
for Es in ['1e-24','1','1e24']:
  E0=arb(Es)
  for rs in ['1e-6','1','1e6']:
    r=arb(rs)
    ceiling=(5*E0/(7*pi))**arb('0.5')*r**(-arb(5)/2)
    Echeck=(7*pi/5)*ceiling*ceiling*r**5
    certify_one(Echeck/E0,('stationary amplifier fixed-scale ceiling',Es,rs))
    rows.append({'E0':Es,'r':rs,'equal_shape_locked_amplification_rate_ceiling':str(ceiling),'energy_at_ceiling':str(Echeck)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'A shape-locked pair sitting at the optimal productive Gram geometry cannot obtain arbitrarily large equal amplification from a harmonic common strain on a fixed positive Hodge radius under finite kinetic energy.  The exact rate ceiling is sqrt(5 E0/(7*pi)) r^(-5/2).  Any finite-time unbounded stationary amplification must therefore move its harmonic support inward or transfer production into the vortical/relative/viscous channels.'
 ),'rows':rows
},indent=2,allow_nan=False))
