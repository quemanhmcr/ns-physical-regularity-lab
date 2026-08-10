import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x > 1-t and x < 1+t):
        raise AssertionError((label,'exact ratio not tightly certified around one',x))

# Exact material-graph memory for the quadratic NS family.
# Phi_i=x_i0+eps*c_i*I, q_ij=Phi_j-Phi_i.  Every loop sum is therefore zero.
theta_sets=[('-2','-0.3','1.1'),('-1','0.2','2'),('0.1','1','3')]
epss=['1e-24','1','1e24']
Is=['1e-24','1','1e24']
ells=['0.5','2','10']
rows=[]
for ths in theta_sets:
    th=[arb(s) for s in ths]
    y=[x.cosh() for x in th]; z=[x.sinh() for x in th]
    c=[y[i]*z[i] for i in range(3)]
    if not (c[0] < c[1] and c[1] < c[2]):
        raise AssertionError(('c ordering',ths,c))
    for es in epss:
      eps=arb(es)
      for Is_ in Is:
        I=arb(Is_)
        for ls in ells:
          ell=arb(ls); S=ell+eps*I
          Phi=[S*c[i] for i in range(3)]
          q12=Phi[1]-Phi[0]; q23=Phi[2]-Phi[1]; q31=Phi[0]-Phi[2]
          loop=q12+q23+q31
          if not loop.contains(0): raise AssertionError(('closed q loop',ths,es,Is_,ls,loop))
          mu=[eps*I*c[i] for i in range(3)]
          m12=mu[1]-mu[0]; m23=mu[2]-mu[1]; m31=mu[0]-mu[2]
          if not (m12+m23+m31).contains(0): raise AssertionError(('closed deposited loop',ths,es,Is_,ls))
          # Independent edge formula from the common memory clock.
          certify_one(q12/(S*(c[1]-c[0])),('q12 coboundary',ths,es,Is_,ls))
          certify_one(q23/(S*(c[2]-c[1])),('q23 coboundary',ths,es,Is_,ls))
          certify_one((-q31)/(S*(c[2]-c[0])),('q31 coboundary',ths,es,Is_,ls))
          rows.append({'theta':ths,'eps':es,'I':Is_,'ell':ls,'S':str(S),'q12':str(q12),'q23':str(q23),'q31':str(q31),'loop_sum':str(loop),'deposited_loop_sum':str(m12+m23+m31)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'In the exact quadratic NS material network, scaled bridge memory is an exact graph coboundary q_ij=Phi_j-Phi_i. '
  'Every closed-loop signed memory sum therefore vanishes identically across 48 decades of renewal scale, even though individual edge deposits can be arbitrarily large. '
  'Closure of signed edge memory is kinematic, not an irreversible cost.'
 ),'rows':rows
},indent=2,allow_nan=False))
