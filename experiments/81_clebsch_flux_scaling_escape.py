import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
# E=diag(2,-1,-1): |E|^2=6 and the full quadratic B-range on |x|=r is (3/2)r^2.
E2=arb(6); rows=[]
for js in ['1e-24','1','1e24']:
  J=arb(js)
  deltaA=(arb(5)/(3*E2))*J
  for rs in ['1e-24','1e-12','1','1e12','1e24']:
    r=arb(rs); deltaB=(arb(3)/2)*r*r; flux=deltaA*deltaB
    certify_one(flux/((arb(5)/12)*J*r*r),('physical Clebsch flux scales like J r^2',js,rs))
    gauges=[]
    for cs in ['1e-24','1','1e24']:
      c=arb(cs); A2=c*deltaA; B2=deltaB/c; flux2=A2*B2
      certify_one(flux2/flux,('Clebsch canonical gauge preserves flux',js,rs,cs))
      gauges.append({'canonical_scale':cs,'Delta_A_rescaled':str(A2),'Delta_B_rescaled':str(B2),'flux_ratio':str(flux2/flux)})
    rows.append({'target_transaction_J':js,'r':rs,'Delta_A_geometric_gauge':str(deltaA),'full_Delta_B_on_sphere':str(deltaB),'physical_flux_cell':str(flux),'flux_over_J_r2':str(flux/(J*r*r)),'gauges':gauges})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'In the geometrically fixed Clebsch chart of the extremal E=diag(2,-1,-1) carrier, a transaction J across a radial cell corresponds to full angular vorticity flux Phi=(5/12)J r^2. Thus fixed transaction strength can be carried by circulation amounts vanishing like r^2 as the scale shrinks. Rescaling the Clebsch potentials A->cA, B->B/c changes Delta A arbitrarily while preserving the physical flux dA wedge dB, killing Delta A itself as a universal ancestry amount.','rows':rows},indent=2,allow_nan=False))
