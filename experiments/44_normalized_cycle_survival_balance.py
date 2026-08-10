import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x > 1-t and x < 1+t):
        raise AssertionError((label,'ratio not tightly certified around one',x))

# Exact normalized survival law for the same quadratic NS family.
# x=e^{-N}; H=c+(L-c)x^2; R2=H^2+4+x^2.
cs=['1e-6','0.01','0.1','1','10','1e6']
Ls=['0.5','2','10']
Ns=['0','1e-12','0.01','0.1','1','2','5','10','30','100','200']
rows=[]
for cs_ in cs:
  c=arb(cs_)
  for Ls_ in Ls:
    L=arb(Ls_)
    for Ns_ in Ns:
      N=arb(Ns_); x=(-N).exp(); x2=x*x
      H=c+(L-c)*x2; R2=H*H+4+x2
      # Rates are per strain time N=g t.
      lambda_D=2*c/H
      sigma_a=arb(1)
      sigma_b=1/(1+x2)
      sigma_R=1 - x2*(2*H*(L-c)+1)/R2
      Tlog_balance=lambda_D-sigma_a-sigma_b-sigma_R
      deficit=sigma_a+sigma_b+sigma_R-lambda_D
      if not (Tlog_balance+deficit).contains(0):
          raise AssertionError(('deficit sign identity',cs_,Ls_,Ns_,Tlog_balance,deficit))

      # Independent direct logarithmic derivative from |T|=x H/[(1+x^2)^1/2 R2^1/2].
      Hlog=-2*(L-c)*x2/H
      R2log=(-4*H*(L-c)*x2-2*x2)/R2
      Tlog_direct=-1 + Hlog + x2/(1+x2) - R2log/2
      if Tlog_direct.contains(0):
          if not Tlog_balance.contains(0): raise AssertionError(('zero direct Tlog mismatch',cs_,Ls_,Ns_))
          ratio='zero/zero'
      else:
          rr=Tlog_balance/Tlog_direct
          certify_one(rr,('normalized survival rate',cs_,Ls_,Ns_))
          ratio=str(rr)

      rows.append({
        'c':cs_,'L':Ls_,'N':Ns_,
        'lambda_D_over_g':str(lambda_D),
        'sigma_a_over_g':str(sigma_a),'sigma_b_over_g':str(sigma_b),'sigma_R_over_g':str(sigma_R),
        'cycle_deficit_over_g':str(deficit),
        'dlog_absT_dN_balance':str(Tlog_balance),'dlog_absT_dN_direct':str(Tlog_direct),
        'balance_over_direct':ratio,
      })

# Large-N exact asymptotic gates: lambda_D->2, endpoint gains->1, bridge rate->1, deficit->1, Tlog->-1.
limits=[]
N=arb(200); x=(-N).exp(); x2=x*x
for cs_ in cs:
  c=arb(cs_)
  for Ls_ in Ls:
    L=arb(Ls_); H=c+(L-c)*x2; R2=H*H+4+x2
    lam=2*c/H; sb=1/(1+x2); sr=1-x2*(2*H*(L-c)+1)/R2
    deficit=1+sb+sr-lam; tlog=-deficit
    for val,target,label in [(lam,arb(2),'lambda_D'),(sb,arb(1),'sigma_b'),(sr,arb(1),'sigma_R'),(deficit,arb(1),'deficit'),(-tlog,arb(1),'minus_Tlog')]:
      certify_one(val/target,('asymptotic '+label,cs_,Ls_),tol='1e-20')
    limits.append({'c':cs_,'L':Ls_,'lambda_D':str(lam),'sigma_b':str(sb),'sigma_R':str(sr),'deficit':str(deficit),'Tlog':str(tlog)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'asymptotic_cases':len(limits),
 'interpretation':(
   'The normalized material-cycle identity d log|T|/dt=lambda_D-sigma_a-sigma_b-sigma_R is independently certified from the closed T formula. '
   'In the exact quadratic NS family, pair-cell renewal approaches two strain rates, the two endpoint magnitude gains consume two, and bridge elongation consumes one more; therefore the cycle renewal deficit approaches one and |T| loses one full exponent per strain time.'
 ),
 'rows':rows,'asymptotic_limits':limits,
},indent=2,allow_nan=False))
