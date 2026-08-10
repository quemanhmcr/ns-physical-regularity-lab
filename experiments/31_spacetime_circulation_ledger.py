import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

# Exact Lamb-Oseen material-spacetime circulation ledger.
# Fixed-radius circles are material. Inventory inside R is
# G_R(t)=G(1-exp(-q)), q=R^2/(4 nu t).
# Outward viscous ancestry current J_R=-dG_R/dt>0.
# Over a time slab, integrated current is exactly
# Loss_R=G_R(t1)-G_R(t2)=G(exp(-q2)-exp(-q1)).
#
# Numerical-observer warning: when q is large, G_R is exponentially close to G.
# Subtracting two Arb balls for G_R can hide the tiny physical loss through
# interval dependency.  The ancestry ledger therefore stores the exponentially
# small tail/current directly rather than reconstructing it from two saturated
# inventories.  This is the physically natural flux variable and fails closed.

base_Rs=['1e-24','1','1e24']
nus=['1e-24','1','1e24']
Gs=['1e-18','1','1e18']
q_pairs=[('100','10'),('10','1'),('3','0.3'),('1','0.1')]
ratios=[arb('0.25'),arb('0.5'),arb(1),arb(2),arb(4),arb(8)]
rows=[]

def exp_tail(q):
    return (-q).exp()

def stable_loss(G,q1,q2):
    # q1>q2. Factor the difference to avoid subtracting saturated inventories.
    # e^-q2 - e^-q1 = e^-q2 [1-e^{-(q1-q2)}].
    return G*exp_tail(q2)*(1-exp_tail(q1-q2))

def annulus_inventory(G,q_inner,q_outer):
    # G_Rout-G_Rin = G(e^-q_inner-e^-q_outer), stable in the vorticity tail.
    if not (q_outer>q_inner):
        raise AssertionError(('annulus q ordering',q_inner,q_outer))
    return G*exp_tail(q_inner)*(1-exp_tail(q_outer-q_inner))

for R0_s in base_Rs:
  R0=arb(R0_s)
  for nu_s in nus:
    nu=arb(nu_s)
    for G_s in Gs:
      G=arb(G_s)
      for q1_s,q2_s in q_pairs:
        q1=arb(q1_s); q2=arb(q2_s)
        if not (q1>q2>0): raise AssertionError('q pair ordering')
        # Base time is defined by q at R0. Same physical t applies to every nested radius.
        t1=R0*R0/(4*nu*q1)
        t2=R0*R0/(4*nu*q2)
        q_at_t1=[]; q_at_t2=[]; losses=[]
        for rr in ratios:
          R=R0*rr
          qi1=R*R/(4*nu*t1)
          qi2=R*R/(4*nu*t2)
          current_loss=stable_loss(G,qi1,qi2)
          if not (current_loss>0):
            raise AssertionError(('nonpositive outward ancestry loss',R0_s,nu_s,G_s,q1_s,q2_s,str(rr),current_loss))
          # Independent antiderivative form of integral J dt after q=c/t:
          # integral J dt = G integral_{q2}^{q1} exp(-q)dq.
          antiderivative_loss=G*(exp_tail(qi2)-exp_tail(qi1))
          # This comparison is only attempted when Arb resolves the raw tail
          # subtraction. If it does not, the stable factorization is the observer.
          if not antiderivative_loss.contains(0):
            if not (current_loss/antiderivative_loss).contains(1):
              raise AssertionError(('current antiderivative mismatch',R0_s,nu_s,G_s,q1_s,q2_s,str(rr),current_loss,antiderivative_loss))
          q_at_t1.append(qi1); q_at_t2.append(qi2); losses.append(current_loss)

        annulus_checks=[]
        sum_delta=arb(0)
        for i in range(len(ratios)-1):
          A1=annulus_inventory(G,q_at_t1[i],q_at_t1[i+1])
          A2=annulus_inventory(G,q_at_t2[i],q_at_t2[i+1])
          delta=A2-A1
          boundary=losses[i]-losses[i+1]
          residual=delta-boundary
          if not residual.contains(0):
            raise AssertionError(('annular spacetime conservation',R0_s,nu_s,G_s,q1_s,q2_s,i,residual,A1,A2,boundary))
          sum_delta += delta
          annulus_checks.append({'inner_ratio':str(ratios[i]),'outer_ratio':str(ratios[i+1]),
                                 'inventory_t1':str(A1),'inventory_t2':str(A2),
                                 'inventory_change':str(delta),'boundary_current_balance':str(boundary),
                                 'residual':str(residual)})

        # Telescoping: sum of all annular changes equals only endpoint-current balance.
        endpoint=losses[0]-losses[-1]
        if not (sum_delta-endpoint).contains(0):
          raise AssertionError(('spacetime telescope failed',R0_s,nu_s,G_s,q1_s,q2_s,sum_delta,endpoint))
        # Direct stable union inventory between the endpoint radii.
        union1=annulus_inventory(G,q_at_t1[0],q_at_t1[-1])
        union2=annulus_inventory(G,q_at_t2[0],q_at_t2[-1])
        union_delta=union2-union1
        if not (union_delta-sum_delta).contains(0):
          raise AssertionError(('union inventory mismatch',R0_s,nu_s,G_s,q1_s,q2_s,union_delta,sum_delta))
        rows.append({'R0':R0_s,'nu':nu_s,'Gamma':G_s,'q1_base':q1_s,'q2_base':q2_s,
                     't2_over_t1':str(t2/t1),'annuli':annulus_checks,
                     'inner_loss':str(losses[0]),'outer_loss':str(losses[-1]),
                     'sum_annular_inventory_change':str(sum_delta),
                     'endpoint_current_balance':str(endpoint),
                     'union_inventory_change':str(union_delta)})

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'cases':len(rows),
  'nested_radii':len(ratios),
  'observer_note':'Exponentially small ancestry loss is represented by the stable tail/current factorization G*exp(-q2)*(1-exp(-(q1-q2))) rather than by subtracting two nearly saturated material-circulation balls.',
  'interpretation':'Exact Lamb-Oseen material worldtubes obey a telescoping ancestry ledger: integrated viscous circulation current through each radius equals the explicit tail loss, annular inventory change is inner current minus outer current, and summing neighboring annuli cancels every internal spacetime boundary exactly. The stable tail representation avoids a false zero caused by subtracting two circulation inventories exponentially close to Gamma. This is the natural anti-double-counting structure for ancestry renewal; it is a conservation law, not yet a universal cost bound.',
  'rows':rows,
},indent=2,allow_nan=False))
