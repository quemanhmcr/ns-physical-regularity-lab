import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x > 1-t and x < 1+t): raise AssertionError((label,x))

theta_sets=[('-2','-0.3','1.1'),('-1','0.2','2'),('0.1','1','3')]
ells=['0.5','2','10']; deposits=['1e-24','1','1e24','1e60']
rows=[]
for ths in theta_sets:
  th=[arb(s) for s in ths]; c=[t.cosh()*t.sinh() for t in th]
  cbar=sum(c)/3; V=sum((ci-cbar)**2 for ci in c)/3
  if not (V>0): raise AssertionError(('nonzero material variance',ths,V))
  for ls in ells:
    ell=arb(ls)
    for ds in deposits:
      dep=arb(ds); S=ell+dep
      x=[S*ci for ci in c]; xb=sum(x)/3
      O=sum((xi-xb)**2 for xi in x)/3
      Ocl=S*S*V; certify_one(O/Ocl,('current occupancy',ths,ls,ds))
      d=[dep*ci for ci in c]; db=sum(d)/3
      M=sum((di-db)**2 for di in d)/3
      Mcl=dep*dep*V
      if ds=='1e-24' or ds=='1' or ds=='1e24' or ds=='1e60': certify_one(M/Mcl,('deposited occupancy',ths,ls,ds))
      mus=[]
      for i,j in [(0,1),(0,2),(1,2)]: mus.append(d[j]-d[i])
      pair=sum(mu*mu for mu in mus)
      certify_one(pair/(9*M),('pairwise-central-moment identity',ths,ls,ds))
      # Gauge translation of all node memory potentials cannot alter the centered occupancy.
      # Keep the observer in the quotient frame: a common gauge shift is tested at a fixed relative scale,
      # rather than recovering tiny centered memory by subtracting an unrelated enormous parent state.
      gauge=dep*arb('1e12'); dg=[di+gauge for di in d]; dgb=sum(dg)/3
      Mg=sum((di-dgb)**2 for di in dg)/3
      certify_one(Mg/M,('gauge-invariant occupancy',ths,ls,ds))
      rows.append({'theta':ths,'ell':ls,'deposit_clock':ds,'V_c':str(V),'current_axial_moment':str(O),'deposited_axial_memory_moment':str(M),'pairwise_edge_memory_square_sum':str(pair),'pair_sum_over_9M':str(pair/(9*M))})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'Signed closed-loop bridge memory cancels exactly, but the centered material occupancy does not. '
  'For the exact hyperbolic triangle the deposited axial second moment is M=(S-ell)^2 V_c, and the sum of squared pairwise deposited separations equals exactly 9M. '
  'This is a physical packet-shape memory and is invariant under a common translation of the node memory potential.'
 ),'rows':rows
},indent=2,allow_nan=False))
