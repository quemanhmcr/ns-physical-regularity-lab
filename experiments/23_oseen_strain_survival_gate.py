import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

# Exact Lamb-Oseen strain survival test.
# At target distance d, q=d^2/(4 nu t), the shear strain magnitude is
# s(q)=Gamma/(2*pi*d^2) F(q), F(q)=1-(1+q)e^{-q}.
# After one INITIAL strain time Delta=1/s(q0),
# q1=q0/[1+8*pi*q0/(Re_Gamma F(q0))], Re_Gamma=Gamma/nu.
# This depends only on q0 and circulation Reynolds number, not on d or nu separately.

def F(q):
    return 1-(1+q)*(-q).exp()

q0s=['0.3','1','3','10','30']
Res=['1e-12','1e-9','1e-6','1e-3','0.01','0.1','1','10','1e2','1e3','1e4','1e6','1e9','1e12']
rows=[]
for q0s_ in q0s:
    q0=arb(q0s_); f0=F(q0)
    if not (f0>0): raise AssertionError(('nonpositive Oseen strain shape',q0s_,f0))
    prev=None
    for Res_ in Res:
        Re=arb(Res_)
        q1=q0/(1+8*pi*q0/(Re*f0))
        f1=F(q1)
        survival=f1/f0
        if not (arb(0)<survival<arb(1)):
            raise AssertionError(('Oseen strain survival outside (0,1)',q0s_,Res_,q1,survival))
        if prev is not None and not (survival>prev):
            raise AssertionError(('strain survival should increase with circulation Reynolds number',q0s_,Res_,prev,survival))
        prev=survival
        small_pred=Re*Re*f0/(128*pi*pi)
        small_ratio=survival/small_pred
        rows.append({'q0':q0s_,'Re_Gamma':Res_,'q1_after_one_initial_strain_time':str(q1),
                     'strain_survival_s1_over_s0':str(survival),
                     'small_Re_asymptotic_ratio':str(small_ratio)})

# Small-Re asymptotic should resolve to 1 without using a fitted exponent.
small=[]
for q0s_ in ['1','3','10']:
    q0=arb(q0s_); f0=F(q0)
    for Res_ in ['1e-6','1e-9','1e-12']:
        Re=arb(Res_)
        q1=q0/(1+8*pi*q0/(Re*f0))
        survival=F(q1)/f0
        pred=Re*Re*f0/(128*pi*pi)
        ratio=survival/pred
        if Res_=='1e-12' and not (abs(ratio-1) < arb('1e-9')):
            raise AssertionError(('small-Re Oseen survival asymptotic not resolved',q0s_,Res_,ratio))
        small.append({'q0':q0s_,'Re_Gamma':Res_,'survival':str(survival),'survival_over_Re2_prediction':str(ratio)})

# Independent scale-cancellation reconstruction from arbitrary d and nu.
scale=[]
for ds in ['1e-24','1e-8','1','1e8','1e24']:
  d=arb(ds)
  for nus in ['1e-24','1','1e24']:
    nu=arb(nus)
    for Res_ in ['0.1','1','10','1e3']:
      Re=arb(Res_); G=Re*nu; q0=arb(10); f0=F(q0)
      t0=d*d/(4*nu*q0)
      s0=G/(2*pi*d*d)*f0
      dt=1/s0
      t1=t0+dt
      q1_physical=d*d/(4*nu*t1)
      q1_reduced=q0/(1+8*pi*q0/(Re*f0))
      if not (q1_physical/q1_reduced).contains(1):
          raise AssertionError(('Oseen survival lost d,nu scale cancellation',ds,nus,Res_,q1_physical/q1_reduced))
      surv=F(q1_physical)/f0
      scale.append({'d':ds,'nu':nus,'Gamma_over_nu':Res_,'q1_ratio_reduced':str(q1_physical/q1_reduced),
                    'strain_survival':str(surv)})

# Concrete gates at q0=10, useful for autopsy but not declared universal constants.
q0=arb(10); f0=F(q0)
def survival_at(Re):
    q1=q0/(1+8*pi*q0/(Re*f0))
    return F(q1)/f0
low=survival_at(arb('0.1'))
unit=survival_at(arb('1'))
high=survival_at(arb('1e4'))
if not (low < arb('1e-4')):
    raise AssertionError(('Re=0.1 donor survived unexpectedly well',low))
if not (unit < arb('0.01')):
    raise AssertionError(('Re=1 donor survived unexpectedly well',unit))
if not (high > arb('0.9')):
    raise AssertionError(('high-Re donor failed to survive one strain time',high))

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'cases':len(rows),
  'small_Re_certificate':small,
  'scale_cancellation_cases':scale,
  'q0_10_gate_samples':{'Re_0.1':str(low),'Re_1':str(unit),'Re_1e4':str(high)},
  'interpretation':'In the exact Lamb-Oseen NS solution, the persistence of strain over one initial strain time depends only on the circulation Reynolds number Gamma/nu and the initial similarity coordinate q0. As Gamma/nu->0 the surviving strain fraction is asymptotic to (Gamma/nu)^2 F(q0)/(128*pi^2), so a vanishing-circulation donor diffuses away before one strain e-fold in this exact model. High-Re circulation survives. This supports, but does not yet prove for arbitrary 3D self-regenerating geometry, a physical survival gate complementary to the closure-collar cost.',
},indent=2,allow_nan=False))
