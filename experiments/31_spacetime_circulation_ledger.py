import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

# Exact Lamb-Oseen material-spacetime circulation ledger.
# Fixed-radius circles are material.  Inventory inside R is
# G_R(t)=G(1-exp(-q)), q=R^2/(4 nu t).
# Outward viscous ancestry current J_R=-dG_R/dt>0.
# Over a time slab, integrated current is exactly
# Loss_R=G_R(t1)-G_R(t2)=G(exp(-q2)-exp(-q1)).
# Nested annular inventories obey Delta A_i = Loss_inner-Loss_outer,
# so all internal spacetime currents telescope exactly.

base_Rs=['1e-24','1','1e24']
nus=['1e-24','1','1e24']
Gs=['1e-18','1','1e18']
q_pairs=[('100','10'),('10','1'),('3','0.3'),('1','0.1')]
ratios=[arb('0.25'),arb('0.5'),arb(1),arb(2),arb(4),arb(8)]
rows=[]

def inv(G,q):
    return G*(1-(-q).exp())

def loss_from_current_integral(G,q1,q2):
    # q decreases as t increases; integrate J dt exactly.
    return G*((-q2).exp()-(-q1).exp())

for R0_s in base_Rs:
  R0=arb(R0_s)
  for nu_s in nus:
    nu=arb(nu_s)
    for G_s in Gs:
      G=arb(G_s)
      for q1_s,q2_s in q_pairs:
        q1=arb(q1_s); q2=arb(q2_s)
        if not (q1>q2>0): raise AssertionError('q pair ordering')
        # Base time is defined by q at R0.  Same physical t applies to every nested radius.
        t1=R0*R0/(4*nu*q1)
        t2=R0*R0/(4*nu*q2)
        inventories1=[]; inventories2=[]; losses=[]
        for rr in ratios:
          R=R0*rr
          qi1=R*R/(4*nu*t1)
          qi2=R*R/(4*nu*t2)
          Gi1=inv(G,qi1); Gi2=inv(G,qi2)
          direct_loss=Gi1-Gi2
          current_loss=loss_from_current_integral(G,qi1,qi2)
          if direct_loss.contains(0):
            raise AssertionError(('unexpected zero loss',R0_s,nu_s,G_s,q1_s,q2_s,str(rr)))
          if not (direct_loss/current_loss).contains(1):
            raise AssertionError(('worldtube current integral mismatch',R0_s,nu_s,G_s,q1_s,q2_s,str(rr),direct_loss/current_loss))
          inventories1.append(Gi1); inventories2.append(Gi2); losses.append(direct_loss)

        annulus_checks=[]
        sum_delta=arb(0)
        for i in range(len(ratios)-1):
          A1=inventories1[i+1]-inventories1[i]
          A2=inventories2[i+1]-inventories2[i]
          delta=A2-A1
          boundary=losses[i]-losses[i+1]
          # The annular inventory can increase or decrease; use an exact residual.
          residual=delta-boundary
          if not residual.contains(0):
            raise AssertionError(('annular spacetime conservation',R0_s,nu_s,G_s,q1_s,q2_s,i,residual))
          sum_delta += delta
          annulus_checks.append({'inner_ratio':str(ratios[i]),'outer_ratio':str(ratios[i+1]),
                                 'inventory_change':str(delta),'boundary_current_balance':str(boundary),
                                 'residual':str(residual)})

        # Telescoping: sum of all annular changes equals only endpoint-current balance.
        endpoint=losses[0]-losses[-1]
        if not (sum_delta-endpoint).contains(0):
          raise AssertionError(('spacetime telescope failed',R0_s,nu_s,G_s,q1_s,q2_s,sum_delta,endpoint))
        # Also reconstruct directly from the union inventory between the two endpoint radii.
        union1=inventories1[-1]-inventories1[0]
        union2=inventories2[-1]-inventories2[0]
        union_delta=union2-union1
        if not (union_delta-sum_delta).contains(0):
          raise AssertionError(('union inventory mismatch',R0_s,nu_s,G_s,q1_s,q2_s,union_delta,sum_delta))
        rows.append({'R0':R0_s,'nu':nu_s,'Gamma':G_s,'q1_base':q1_s,'q2_base':q2_s,
                     't2_over_t1':str(t2/t1),'annuli':annulus_checks,
                     'sum_annular_inventory_change':str(sum_delta),
                     'endpoint_current_balance':str(endpoint),
                     'union_inventory_change':str(union_delta)})

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'cases':len(rows),
  'nested_radii':len(ratios),
  'interpretation':'Exact Lamb-Oseen material worldtubes obey a telescoping ancestry ledger: integrated viscous circulation current through each radius equals the explicit loss of material circulation inside it, annular inventory change is inner current minus outer current, and summing neighboring annuli cancels every internal spacetime boundary exactly. This is the natural anti-double-counting structure for ancestry renewal; it is a conservation law, not yet a universal cost bound.',
  'rows':rows,
},indent=2,allow_nan=False))
