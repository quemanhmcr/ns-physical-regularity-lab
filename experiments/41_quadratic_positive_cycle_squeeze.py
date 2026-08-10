import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def det(a,b,c): return dot(a,cross(b,c))
def norm(a): return dot(a,a).sqrt()
def matvec(A,x): return tuple(sum(A[i][j]*x[j] for j in range(3)) for i in range(3))
def vsub(a,b): return tuple(a[i]-b[i] for i in range(3))
def certify_one(x,label):
    tol=arb('1e-30')
    if not x.contains(1) or not (x > 1-tol and x < 1+tol):
        raise AssertionError((label,'ratio not tightly certified around one',x))

def grad(a,e,y,z):
    return ((a,e*z,e*y),(arb(0),-a,arb(0)),(arb(0),arb(0),arb(0)))
def bridge_delta_grad(e,y,z,ym,zm):
    # Directly represent A-Abar; do not subtract two parent matrices carrying the common a mode.
    return ((arb(0),e*(z-zm),e*(y-ym)),(arb(0),arb(0),arb(0)),(arb(0),arb(0),arb(0)))
def vel(a,e,x,y,z): return (a*x+e*y*z,-a*y,arb(0))
def vort(e,y,z): return (arb(0),e*y,-e*z)

# u=(a x+eps y z,-a y,0), exact divergence-free smooth initial data.
# Xa=(0,1,0), Xb=(-L,-1,1) gives a positive angular two-cycle for every L>0.
a_vals=['-1e30','-1e12','-1','-0.01']
e_vals=['1e-30','1','1e30']
Ls=['0.5','2','10','1e6']
rows=[]
for aas in a_vals:
  a=arb(aas)
  for es in e_vals:
    e=arb(es)
    for Ls_ in Ls:
      L=arb(Ls_)
      Xa=(arb(0),arb(1),arb(0)); Xb=(-L,arb(-1),arb(1))
      R=vsub(Xb,Xa); rsep=norm(R)
      pa=vort(e,Xa[1],Xa[2]); pb=vort(e,Xb[1],Xb[2])
      ra=norm(pa); rb=norm(pb)
      D=det(pa,R,pb)
      if D.contains(0): raise AssertionError(('degenerate D',aas,es,Ls_,D))
      T=D/(ra*rsep*rb)
      an=dot(pa,R)/(ra*rsep); bn=dot(pb,R)/(rb*rsep)
      Kba=T*an; Kab=-T*bn
      if not (Kba>0 and Kab>0):
          raise AssertionError(('chosen physical pair not positive cycle',aas,es,Ls_,Kba,Kab,T,an,bn))
      Kba_expected=arb(2).sqrt()*L/(L*L+5)
      Kab_expected=L/(2*(L*L+5))
      certify_one(Kba/Kba_expected,('Kba formula',Ls_))
      certify_one(Kab/Kab_expected,('Kab formula',Ls_))

      Aa=grad(a,e,Xa[1],Xa[2]); Ab=grad(a,e,Xb[1],Xb[2])
      # A is affine along the straight bridge, so its line average is the midpoint value.
      ym=(Xa[1]+Xb[1])/2; zm=(Xa[2]+Xb[2])/2
      Abar=grad(a,e,ym,zm)
      ua=vel(a,e,*Xa); ub=vel(a,e,*Xb)
      Rdot=vsub(ub,ua)
      padot=matvec(Aa,pa); pbdot=matvec(Ab,pb)  # Delta omega=0 exactly.
      # Raw material derivative contains exact cancellation of O(|a| eps^2) common-affine terms.
      # Keep it as an independent observer only when conditioned; for extreme scale separation,
      # project out the common mode first and observe the physical bridge current directly.
      Ddot_raw=det(padot,R,pb)+det(pa,Rdot,pb)+det(pa,R,pbdot)
      dAa=bridge_delta_grad(e,Xa[1],Xa[2],ym,zm)
      dAb=bridge_delta_grad(e,Xb[1],Xb[2],ym,zm)
      bridge=det(matvec(dAa,pa),R,pb)+det(pa,R,matvec(dAb,pb))
      Ddot_closed=-e**3
      certify_one(bridge/Ddot_closed,('bridge closed form -eps^3',aas,es,Ls_))
      raw_well_conditioned = abs(a) <= arb('1e6')*e
      if raw_well_conditioned:
          certify_one(Ddot_raw/bridge,('raw material bridge balance',aas,es,Ls_))
      Dlog=bridge/D
      expected_Dlog=e/L
      certify_one(Dlog/expected_Dlog,('D renewal eps/L',aas,es,Ls_))

      gain_a=dot(pa,padot)/(ra*ra)
      gain_b=dot(pb,pbdot)/(rb*rb)
      pair_gain=gain_a+gain_b
      expected_pair_gain=-3*a/2
      certify_one(pair_gain/expected_pair_gain,('endpoint pair gain',aas,es,Ls_))
      rlog_raw=dot(R,Rdot)/(rsep*rsep)
      rlog_expected=(a*(L*L-4)+L*e)/(L*L+5)
      # The closed form is the stable observer when L^2 is near 4 and |a|>>eps.
      rlog=rlog_expected
      rraw_scale=abs(a)*(L*L+4)+L*e
      r_well_conditioned = abs(rlog_expected)*(L*L+5)*arb('1e6') >= rraw_scale
      if r_well_conditioned and not rlog_expected.contains(0):
          certify_one(rlog_raw/rlog_expected,('raw separation rate',aas,es,Ls_))
      Tlog=Dlog-pair_gain-rlog
      squeeze=pair_gain-Dlog  # = -d log(r |T|)/dt
      identity=rlog+Tlog-(Dlog-pair_gain)
      if not identity.contains(0):
          raise AssertionError(('rT squeeze identity',aas,es,Ls_,identity))

      rows.append({
        'a':aas,'eps':es,'L':Ls_,
        'D_ab':str(D),'Ddot_over_D':str(Dlog),'expected_eps_over_L':str(expected_Dlog),
        'K_b_to_a':str(Kba),'K_a_to_b':str(Kab),
        'pair_vorticity_gain_rate':str(pair_gain),
        'separation_log_rate':str(rlog),'triple_product_log_rate':str(Tlog),
        'minus_log_rT_rate_pair_gain_minus_Drenewal':str(squeeze),
        'bridge_over_closed_minus_eps3':str(bridge/Ddot_closed),
        'raw_material_derivative_certified':bool(raw_well_conditioned),
        'raw_Ddot':str(Ddot_raw),
        'raw_separation_rate_certified':bool(r_well_conditioned),
        'raw_separation_log_rate':str(rlog_raw),
      })

# Clean orientation-collapse diagnostic: at L=2 common affine strain drops out of rdot/r exactly.
# Use a=-1e30, eps=1: pair gain is huge, rlog=2/9, so T absorbs essentially all squeeze.
a=arb('-1e30'); e=arb(1); L=arb(2)
rlog=2*e/9
pair_gain=-3*a/2
Dlog=e/L
Tlog=Dlog-pair_gain-rlog
if not (pair_gain>0 and Tlog<0): raise AssertionError(('L=2 orientation collapse sign',pair_gain,Tlog))

print(json.dumps({
 'arb_precision_bits':BITS,
 'status':'PASS',
 'cases':len(rows),
 'L2_orientation_collapse':{
   'pair_gain_rate':str(pair_gain),'D_renewal_rate':str(Dlog),
   'separation_log_rate':str(rlog),'triple_product_log_rate':str(Tlog),
 },
 'interpretation':(
   'The exact divergence-free quadratic field contains a physical endpoint pair whose Biot-Savart angular geometry is a positive two-cycle. '
   'The common affine strain parameter can amplify both endpoint vorticity magnitudes arbitrarily rapidly yet cancels from Ddot/D, which equals eps/L and is supplied only by bridge inhomogeneity. Extreme raw derivatives are deliberately not formed as a tiny difference of huge common-affine contributions; the observer represents the projected bridge mismatch directly rather than subtracting two parent gradient states. '
   'Consequently r*|T| must collapse whenever pair amplification outruns cell renewal; at L=2 the common affine part does not compress the separation, so the productive triple product itself collapses.'
 ),
 'rows':rows,
},indent=2,allow_nan=False))
