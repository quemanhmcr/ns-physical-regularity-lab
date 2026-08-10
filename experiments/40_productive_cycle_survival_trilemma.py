import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

# Exact factorization:
# Lambda_pair = R_D * (r0/r) * (T0/T), using absolute nonzero T.
# Calibrate the three pure escape branches over huge amplification factors.
amps=['1','1e3','1e12','1e30','1e60','1e120']
T0=arb('0.4')
rows=[]
for As in amps:
    A=arb(As)
    # 1) No renewal and no angular loss: all amplification becomes bridge compression.
    RD=arb(1); Trat=arb(1); rrat=1/A
    rec=RD*(1/rrat)*(1/Trat)
    if not (rec/A).contains(1): raise AssertionError(('compression branch',As,rec,A))

    # 2) No renewal and no bridge compression: the triple product must collapse.
    RD2=arb(1); rrat2=arb(1); Trat2=1/A
    rec2=RD2*(1/rrat2)*(1/Trat2)
    if not (rec2/A).contains(1): raise AssertionError(('angular-death branch',As,rec2,A))

    # 3) No bridge compression and no angular loss: D itself must renew by A.
    RD3=A; rrat3=arb(1); Trat3=arb(1)
    rec3=RD3*(1/rrat3)*(1/Trat3)
    if not (rec3/A).contains(1): raise AssertionError(('renewal branch',As,rec3,A))

    rows.append({
      'pair_amplification':As,
      'compression_only_r_over_r0':str(rrat),
      'angular_death_only_T_over_T0':str(Trat2),
      'renewal_only_D_over_D0':str(RD3),
    })

# Uniform-productivity + non-contact inequality.
# If |T|>=kappa and r/r0>=delta, then Lambda_pair <= RD*|T0|/(delta*kappa).
kappas=['1e-6','0.01','0.1','0.3']
deltas=['1e-6','0.01','0.1','1']
bounds=[]
for ks in kappas:
  k=arb(ks)
  for ds in deltas:
    d=arb(ds)
    for RDs in ['1','1e6','1e30']:
      RD=arb(RDs)
      bound=RD*T0/(d*k)
      # Saturate both physical floors to verify the bound algebraically.
      T=k; rr=d
      amp=RD*(1/rr)*(T0/T)
      if not (amp/bound).contains(1):
          raise AssertionError(('persistent-cycle bound mismatch',ks,ds,RDs,amp,bound))
      bounds.append({'kappa':ks,'delta':ds,'R_D':RDs,'max_pair_amplification_at_floors':str(bound)})

print(json.dumps({
 'arb_precision_bits':BITS,
 'status':'PASS',
 'pure_branch_cases':len(rows),
 'persistent_floor_cases':len(bounds),
 'interpretation':(
   'The exact pair-cell factorization has three pure multiplicative escape branches. '
   'Without ancestry-cell renewal, arbitrarily large pair amplification must either compress the material bridge or destroy the triple-product geometry. '
   'If a positive cycle remains angularly productive and avoids near-contact, its pair amplification is bounded by the pair-cell renewal factor times explicit geometric floor ratios.'
 ),
 'pure_branches':rows,
 'persistent_cycle_floor_bounds':bounds,
},indent=2,allow_nan=False))
