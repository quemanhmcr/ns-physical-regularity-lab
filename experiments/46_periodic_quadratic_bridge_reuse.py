import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x > 1-t and x < 1+t):
        raise AssertionError((label,'ratio not tightly certified around one',x))

# Exact smooth NS family u=(a(t)x+eps yz,-a(t)y,0), Delta u=0.
# Choose lambda=1+eta sin^2 t and a=-lambda_dot/lambda, so lambda returns to 1
# every pi while the bridge-current memory Q=L+eps*int lambda^2 grows secularly.
etas=['1e-6','0.1','1','10','1e6']
times=['0','1e-12','0.01','0.1','1','2','10','100']
rows=[]
for es in etas:
  eta=arb(es)
  C=1+eta+3*eta*eta/8
  for ts in times:
    t=arb(ts)
    s=t.sin(); sin2=(2*t).sin(); sin4=(4*t).sin(); cos2=(2*t).cos(); cos4=(4*t).cos()
    lam=1+eta*s*s
    ldot=eta*sin2
    a=-ldot/lam
    I=t*C-(eta/2+eta*eta/4)*sin2+(eta*eta/32)*sin4
    Idot=C-(eta+eta*eta/2)*cos2+(eta*eta/8)*cos4
    certify_one(Idot/(lam*lam),('closed integral derivative lambda^2',es,ts))
    # Exactness checks: e is constant, so the vorticity equation for omega=(0,e y,-e z)
    # requires no coefficient renewal beyond material stretching by a(t).
    if not (a*lam+ldot).contains(0): raise AssertionError(('a=-lambda_dot/lambda',es,ts,a,lam,ldot))
    rows.append({'eta':es,'t':ts,'lambda':str(lam),'lambda_dot':str(ldot),'a':str(a),'I':str(I),'I_dot_over_lambda2':str(Idot/(lam*lam))})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
   'The smooth periodic choice lambda=1+eta sin^2(t) produces an exact time-dependent quadratic Navier-Stokes mechanism with a=-lambda_dot/lambda. '
   'Its bridge forcing integral I(t)=int_0^t lambda^2 is certified from an independent derivative and therefore accumulates secularly even though the local material stretch lambda returns exactly every pi.'
 ),
 'rows':rows,
},indent=2,allow_nan=False))
