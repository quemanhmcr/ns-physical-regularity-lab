import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rt2=arb(2).sqrt()

def qnorm2(Q): return sum(Q[i][j]*Q[i][j] for i in range(3) for j in range(3))
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))

# Unit stationary capacity amplifier S_unit=S_*/sigma.
c=arb(3)/(2*rt2)
S=((arb(1),arb(0),c),(arb(0),arb(1),c),(c,c,-arb(2)))
S2=qnorm2(S); certify_one(S2/(arb(21)/2),('capacity amplifier unit norm',))
rows=[]
for qs in ['1e-24','1','1e24']:
  q=arb(qs); Q2=q*q*S2
  floor=(20*pi/9)*Q2
  closed=(70*pi/3)*q*q
  certify_one(floor/closed,('capacity pure-production angular floor',qs))
  for rs in ['1e-6','1','1e6']:
    r=arb(rs)
    Gamma_norm=r*r*(Q2.sqrt())
    rows.append({'production_channel_q':qs,'r':rs,'Q_Frobenius_squared':str(Q2),'sharp_tangential_vorticity_sphere_floor':str(floor),'closed_70pi_over3_q2':str(closed),'transaction_circulation_Frobenius_norm_r2Q':str(Gamma_norm)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'S_unit_Frobenius_squared':str(S2),
 'interpretation':(
  'For the pair-frame ray that produces equal positive magnitude amplification while leaving all three Gram-shape channels zero at the capacity geometry, Q=q S_unit and |S_unit|^2=21/2. The universal sharp transaction inequality therefore specializes to integral|n cross omega|^2 >= (70*pi/3)q^2 on every supplying sphere. The natural circulation-dimensional size is r^2|Q|.'
 ),'rows':rows
},indent=2,allow_nan=False))
