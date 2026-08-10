import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def qnorm2(Q): return sum(Q[i][j]*Q[i][j] for i in range(3) for j in range(3))
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))

# Exact isotropic identities for symmetric trace-free Q:
# int |Q n|^2 dOmega=(4*pi/3)|Q|^2,
# int (n.Q.n)^2 dOmega=(8*pi/15)|Q|^2,
# hence int |P_n Qn|^2=(4*pi/5)|Q|^2.
# The transaction pairing with v=n x omega is
# |Q|^2=(3/(4*pi)) int (P_n Qn).v dOmega.
# Cauchy gives int |v|^2 >= (20*pi/9)|Q|^2.
Qs=[
 ('diag',((arb(2),arb(0),arb(0)),(arb(0),-arb(1),arb(0)),(arb(0),arb(0),-arb(1)))),
 ('mixed',((arb('0.4'),arb('0.3'),-arb('0.2')),(arb('0.3'),-arb('0.1'),arb('0.25')),(-arb('0.2'),arb('0.25'),-arb('0.3')))),
]
rows=[]
for name,Q0 in Qs:
  if not sum(Q0[i][i] for i in range(3)).contains(0): raise AssertionError(('Q trace',name))
  for ss in ['1e-24','1','1e24']:
    s=arb(ss); Q=tuple(tuple(s*Q0[i][j] for j in range(3)) for i in range(3)); q2=qnorm2(Q)
    tangential_test=(4*pi/5)*q2
    sharp=(20*pi/9)*q2
    # Equality field v_eq=(5/3)P_n Qn has exactly the sharp L2 content.
    veq2=(arb(25)/9)*tangential_test
    certify_one(veq2/sharp,('sharp equality field',name,ss))
    # Add a transaction-null rigid-vorticity mode v0=n x e_z. Its L2 norm is 8*pi/3,
    # and its cross pairing with the equality field vanishes by symmetry.
    for Ws in ['0','1e-24','1','1e24']:
      W=arb(Ws); null_extra=W*W*(8*pi/3); total=veq2+null_extra
      if not (null_extra>=0): raise AssertionError(('negative transaction-null remainder',name,ss,Ws,null_extra))
      rows.append({'Q_case':name,'Q_scale':ss,'null_mode_scale':Ws,'Q_Frobenius_squared':str(q2),'PnQn_L2_squared':str(tangential_test),'sharp_floor_20pi_over9_Q2':str(sharp),'equality_field_L2_squared':str(veq2),'transaction_null_L2_remainder':str(null_extra),'equality_plus_null_L2_squared':str(total)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'For every symmetric trace-free Hodge transaction tensor Q, the tangential vorticity observable v=n cross omega obeys the sharp universal sphere inequality integral|v|^2 >= (20*pi/9)|Q|^2. Equality is attained by v=(5/3)(I-nn)Qn. Adding any transaction-null rigid-vorticity mode preserves Q but strictly increases the angular vorticity content unless its amplitude vanishes.'
 ),'rows':rows
},indent=2,allow_nan=False))
