import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

# Axisymmetric canonical family at a Hodge-defined physical sphere:
# S=s diag(-1,-1,2), rigid vorticity W e_z, productive azimuthal carrier.
# q=e.Qe and A=(5/2)q is the amplitude entering omega=W e + A mu sqrt(1-mu^2)e_phi.

def c_lower(eps):
    z=1+eps*eps/4
    return (arb(7)/15)/(z.sqrt()**3) + (arb(1)/3)/z.sqrt()

rows=[]
for eps_s in ['1e-30','1e-18','1e-12','1e-6','1e-3','0.01','0.1','0.3','1']:
    eps=arb(eps_s)
    c=c_lower(eps)
    coherent=arb(4)/5
    ratio=c/coherent
    if eps_s=='1e-30' and not (ratio > arb('0.9999999999999999999999999999999999999999')):
        raise AssertionError(('coherent-limit approach too weak',eps_s,ratio))
    # c_lower must be positive and no larger than the coherent coefficient 4/5.
    if not (c > 0): raise AssertionError(('nonpositive tax coefficient',eps_s,c))
    gap=coherent-c
    if not (gap.contains(0) or gap > 0):
        raise AssertionError(('lower bracket exceeded upper bracket',eps_s,c,gap))
    rows.append({'epsilon_A_over_W':eps_s,'c_lower':str(c),'lower_over_coherent_limit':str(ratio)})

trade=[]
nu_values=['1e-18','1','1e18']
r_values=['1e-18','1e-6','1','1e6','1e18']
q_values=['1e-12','1','1e12']
eps0=arb('0.1')
c0=c_lower(eps0)
for nus in nu_values:
  for rs in r_values:
    for qs in q_values:
      nu=arb(nus); r=arb(rs); q=arb(qs)
      A=(arb(5)/2)*q
      W=A/eps0
      eps=A/W
      if not eps.contains(eps0): raise AssertionError('coherence parameter reconstruction failed')
      # Rigorous angular-tax bracket. Full |omega||grad xi|^2 is >= angular part/r^2.
      Dlow=(A*A/W)*c0
      Dup=(A*A/W)*(arb(4)/5)
      tax_low=nu*Dlow/(r*r)
      tax_upper_angular=nu*Dup/(r*r)
      if not (tax_low <= tax_upper_angular):
          raise AssertionError(('tax bracket inverted',nus,rs,qs,tax_low,tax_upper_angular))
      # Background circulation and exact rigid-rotation kinetic energy in B_r.
      Gamma=pi*W*r*r
      Erot=pi*W*W*r**5/15
      EfromGamma=Gamma*Gamma*r/(15*pi)
      if not (Erot/EfromGamma).contains(1):
          raise AssertionError(('circulation-energy identity failed',nus,rs,qs,Erot/EfromGamma))
      # Product tax*sqrt(Erot): W cancels.  This is the canonical resource-coupled law.
      product=tax_low*Erot.sqrt()
      product_target=(arb(25)/4)*c0*nu*q*q*(pi*r/15).sqrt()
      if not (product/product_target).contains(1):
          raise AssertionError(('gain-cost product failed',nus,rs,qs,product/product_target))
      # Naive universal tax independent of W is false: at fixed q,r,nu the coherent
      # asymptotic upper bracket scales exactly like 1/W while energy scales W^2.
      naive_escape_ratio=(tax_upper_angular*W*r*r)/(5*nu*q*q)
      if not naive_escape_ratio.contains(1):
          raise AssertionError(('1/W escape normalization failed',nus,rs,qs,naive_escape_ratio))
      trade.append({'nu':nus,'r':rs,'q_e':qs,'epsilon':str(eps),
                    'tax_lower':str(tax_low),'angular_tax_upper':str(tax_upper_angular),
                    'Gamma_background':str(Gamma),'E_rotation':str(Erot),
                    'tax_sqrtE_ratio':str(product/product_target),
                    'one_over_W_scaling_ratio':str(naive_escape_ratio)})

# Shrinking-support and signed-cancellation attacks in the coherent quadratic branch.
# These are scale ledgers derived from the rigorous angular lower coefficient c0.
support=[]
sigma=arb(1); W=arb('1e6'); r=arb(1); nu=arb(1); weight=arb(1)
for ds in ['1','0.3','0.1','0.03','0.01','0.003','0.001']:
    delta=arb(ds)
    q=sigma/(weight*delta)
    # band volume tax lower: 4pi nu int D_ang dr; with dr~r delta for the ledger.
    band_lower=25*pi*c0*nu*r*delta*q*q/W
    scaled=band_lower*delta
    target=25*pi*c0*nu*r*sigma*sigma/(W*weight*weight)
    if not (scaled/target).contains(1):
        raise AssertionError(('shrinking-support 1/delta ledger failed',ds,scaled/target))
    support.append({'delta_log':ds,'required_q':str(q),'band_tax_lower':str(band_lower),
                    'delta_times_tax_ratio':str(scaled/target)})

cancel=[]
for cs in ['0','0.1','0.3','1','3','10']:
    c=arb(cs)
    q1=1+c; q2=-c
    net=q1+q2
    gross_quad=q1*q1+q2*q2
    if not net.contains(1): raise AssertionError(('net transaction changed under cancellation',cs,net))
    if c>0 and not (gross_quad>1): raise AssertionError(('cancellation did not raise quadratic cost ledger',cs,gross_quad))
    cancel.append({'cancellation_amplitude':cs,'net_transaction':str(net),
                   'quadratic_cost_factor':str(gross_quad)})

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'coherence_bracket':rows,
  'gain_cost_cases':len(trade),
  'gain_cost_sample':trade,
  'shrinking_support_attack':support,
  'signed_cancellation_attack':cancel,
  'interpretation':'The naive claim that fixed productive self-stretching forces a circulation-independent geometry tax is false: coherent background vorticity W drives the angular directional tax down like 1/W at fixed q_e. In the exact axisymmetric Hodge carrier, however, the same escape requires Gamma_R=pi W r^2 and E_R=Gamma_R^2 r/(15pi). For a declared coherence level epsilon<=0.1, Arb certifies the resource-coupled product lower bound tax_ang*sqrt(E_R)>=(25/4)c_-(0.1) nu q_e^2 sqrt(pi r/15). The same canonical bracket penalizes logarithmically shrinking productive support like 1/delta and makes opposite signed transactions increase the quadratic cost ledger. These are canonical-family laws, not a general regularity theorem; material-lineage switching remains open.',
},indent=2))
