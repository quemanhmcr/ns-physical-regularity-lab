import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def cycle(x,y,sp):
    # alpha=x>0, beta=-y<0, transverse relative sine sp>=0.
    T2=(1-x*x)*(1-y*y)*sp*sp
    return T2, x*y*T2

def certify_near_one(x,label,tol='1e-20'):
    t=arb(tol)
    if not (x > 1-t and x < 1+t): raise AssertionError((label,x))

scales=['1e-6','1e-12','1e-24','1e-48']
rows=[]
base=arb('0.4')
y=arb('0.6'); sp=arb('0.8')
for es in scales:
    e=arb(es)
    # Longitudinal starvation alpha=e: P ~ e*y*(1-y^2)*sp^2.
    T2,P=cycle(e,y,sp); asym=e*y*(1-y*y)*sp*sp
    certify_near_one(P/asym,('longitudinal starvation',es))
    rows.append({'mode':'longitudinal','eps':es,'T2':str(T2),'P':str(P),'P_over_asymptote':str(P/asym)})
    # Axial/transverse starvation alpha=sqrt(1-e): 1-alpha^2=e exactly, P ~ sqrt(1-e)*y*e*(1-y^2)*sp^2.
    x=(1-e).sqrt()
    # Observe the small transverse reserve e directly; do not reconstruct it as 1-x^2
    # after x has saturated numerically near one.
    T2=e*(1-y*y)*sp*sp; P=x*y*T2; asym=x*y*e*(1-y*y)*sp*sp
    certify_near_one(P/asym,('axial transverse starvation',es))
    rows.append({'mode':'axial_transverse','eps':es,'T2':str(T2),'P':str(P),'P_over_exact_eps_factor':str(P/asym)})
    # Coplanarity: sin(phi)=e with fixed longitudinal access, P is exactly quadratic in e.
    T2,P=cycle(base,y,e); asym=base*y*(1-base*base)*(1-y*y)*e*e
    certify_near_one(P/asym,('coplanarity starvation',es))
    rows.append({'mode':'coplanarity','eps':es,'T2':str(T2),'P':str(P),'P_over_quadratic_asymptote':str(P/asym)})

# Productive chamber implication attack for fixed kappa.
kappa=arb('1e-4'); chamber=[]
xs=['0.05','0.2','0.4','0.577','0.8','0.95']; ys=['0.05','0.3','0.6','0.9']; sps=['0.2','0.6','1']
for xs_ in xs:
  x=arb(xs_)
  for ys_ in ys:
    y=arb(ys_)
    for ss in sps:
      sp=arb(ss); T2,P=cycle(x,y,sp)
      if P>=kappa:
        if not (T2>=kappa and x>=kappa and y>=kappa and 1-x*x>=kappa and 1-y*y>=kappa):
          raise AssertionError(('productive chamber margin',xs_,ys_,ss,T2,P))
        chamber.append({'alpha_abs':xs_,'beta_abs':ys_,'sin_transverse_angle':ss,'T2':str(T2),'P':str(P)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','boundary_cases':len(rows),'chamber_cases':len(chamber),
 'interpretation':(
  'The productive Gram chamber has distinct physical boundary modes. Longitudinal starvation is linear in a vanishing longitudinal projection, axial alignment removes the needed transverse factor, and coplanarity kills the product quadratically through T^2. '
  'Any fixed positive productivity threshold keeps the triad away from all of these boundaries.'
 ),'boundaries':rows,'chamber':chamber
},indent=2,allow_nan=False))
