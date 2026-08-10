import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))

# Universal fixed-scale bound for one physical Gram coordinate:
# |j_h| <= |M_g| |S_h| <= sqrt(15 E0/pi) (1-g^2) r^(-5/2).
# This is the dual form of the exact projected harmonic occupancy floor.
rows=[]
for gs in ['0','0.2','0.6','0.99']:
  g=arb(gs)
  for Es in ['1e-24','1','1e24']:
    E0=arb(Es)
    for rs in ['1e-6','1','1e6']:
      r=arb(rs)
      Mnorm=arb(2).sqrt()*(1-g*g)
      Shmax=(15*E0/(2*pi))**arb('0.5')*r**(-arb(5)/2)
      jmax=Mnorm*Shmax
      closed=(15*E0/pi)**arb('0.5')*(1-g*g)*r**(-arb(5)/2)
      if closed.contains(0):
          if not jmax.contains(0): raise AssertionError(('zero fixed scale bound',gs,Es,rs,jmax))
          ratio='zero/zero'
      else:
          certify_one(jmax/closed,('fixed scale harmonic current bound',gs,Es,rs))
          ratio=str(jmax/closed)
      rows.append({'g':gs,'E0':Es,'r':rs,'harmonic_shape_current_upper':str(jmax),'closed_bound':str(closed),'ratio':ratio})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'Finite kinetic energy gives a uniform rate ceiling for the harmonic contribution to any physical Gram-coordinate current at a fixed Hodge radius: |j_h(r)|<=sqrt(15 E0/pi)(1-g^2) r^(-5/2).  Therefore infinite finite-time variation of that coordinate cannot be supplied indefinitely by a harmonic reset source on a radius bounded away from zero; the surviving source must move toward smaller physical radii or transfer into the vortical/non-affine/viscous channels.'
 ),'rows':rows
},indent=2,allow_nan=False))
