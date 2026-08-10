import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

# Angular part of the pairwise Biot-Savart stretching kernel.
# Target direction a, donor direction b, separation n from a to b:
#   K_ba=(a.n) a.(n x b)
# Reverse transaction:
#   K_ab=(b.(-n)) b.((-n) x a)=-(b.n) a.(n x b).
# A common positive strength/distance prefactor is omitted; the experiment
# isolates the exact signed angular geometry of the directed transaction graph.

sqrt2=arb(2).sqrt()
sqrt3=arb(3).sqrt()

# One-way exact configuration.
# n=e_z, a=(1/sqrt2,0,1/sqrt2), b=(0,-1,0).
T_one=1/sqrt2
an_one=1/sqrt2
bn_one=arb(0)
K_ba_one=T_one*an_one
K_ab_one=-T_one*bn_one
if not K_ba_one.contains(arb('0.5')):
    raise AssertionError(('one-way forward edge not 1/2',K_ba_one))
if not K_ab_one.contains(0):
    raise AssertionError(('one-way reverse edge not zero',K_ab_one))

# Symmetric positive two-cycle family parameterized by c=cos(theta).
cs=['1e-12','0.01','0.1','0.25','0.5','0.577350269189625764509148780501957456','0.75','0.9','0.99']
rows=[]
for cs_ in cs:
    c=arb(cs_)
    s2=1-c*c
    if not (s2>0):
        raise AssertionError(('invalid symmetric family cosine',cs_,s2))
    T=s2
    an=c
    bn=-c
    Kba=T*an
    Kab=-T*bn
    if not (Kba>0 and Kab>0):
        raise AssertionError(('mutual-stretching pair lost positive sign',cs_,Kba,Kab))
    if not (Kba/Kab).contains(1):
        raise AssertionError(('symmetric two-cycle lost equal edge weights',cs_,Kba,Kab))
    product=Kba*Kab
    product_identity=product/(-T*T*an*bn)
    if not product_identity.contains(1):
        raise AssertionError(('pair-product identity failed',cs_,product_identity))
    rows.append({
        'cos_theta':cs_,
        'triple_product_T':str(T),
        'a_dot_n':str(an),
        'b_dot_n':str(bn),
        'K_b_to_a':str(Kba),
        'K_a_to_b':str(Kab),
        'edge_product':str(product),
        'product_identity_ratio':str(product_identity),
    })

# Exact maximum f(c)=c(1-c^2) occurs at c=1/sqrt(3).
cstar=1/sqrt3
fstar=cstar*(1-cstar*cstar)
expected=2/(3*sqrt3)
if not (fstar/expected).contains(1):
    raise AssertionError(('mutual-stretch maximum identity failed',fstar,expected))
# Exact derivative f'=1-3c^2 vanishes there.
fprime=1-3*cstar*cstar
if not fprime.contains(0):
    raise AssertionError(('mutual-stretch maximum derivative excludes zero',fprime))
# Compare against the sampled shoulders without relying on a numerical optimizer.
for cs_ in ['0.25','0.5','0.75','0.9']:
    c=arb(cs_); f=c*(1-c*c)
    if not (f<fstar):
        raise AssertionError(('symmetric maximum shoulder attack failed',cs_,f,fstar))

# Sign taxonomy from the exact product formula.  Positive same-sign cycles need
# opposite longitudinal projections.  Same longitudinal signs force opposite edge signs.
sign_cases=[]
for T_s in ['-2','-0.1','0.1','2']:
  T=arb(T_s)
  for an_s,bn_s in [('0.5','-0.25'),('-0.5','0.25'),('0.5','0.25'),('-0.5','-0.25')]:
    an=arb(an_s); bn=arb(bn_s)
    Kba=T*an; Kab=-T*bn
    prod=Kba*Kab
    expected_prod=-T*T*an*bn
    if not (prod/expected_prod).contains(1):
        raise AssertionError(('sign taxonomy product identity',T_s,an_s,bn_s))
    opposite_longitudinal=(an*bn<0)
    same_edge_sign=(prod>0)
    if bool(opposite_longitudinal) != bool(same_edge_sign):
        raise AssertionError(('same-sign cycle iff opposite longitudinal tilt failed',T_s,an_s,bn_s,prod))
    sign_cases.append({
        'T':T_s,'a_dot_n':an_s,'b_dot_n':bn_s,
        'K_b_to_a':str(Kba),'K_a_to_b':str(Kab),
        'same_edge_sign':bool(same_edge_sign),
        'opposite_longitudinal_projections':bool(opposite_longitudinal),
    })

# Coplanar/null control: T=0 kills both directions regardless of longitudinal tilt.
for an_s,bn_s in [('1','-1'),('1','1'),('1e30','-1e30')]:
    an=arb(an_s); bn=arb(bn_s); T=arb(0)
    if not ((T*an).contains(0) and (-T*bn).contains(0)):
        raise AssertionError('coplanar triple-product null failed')

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'one_way_edge':{
      'K_b_to_a':str(K_ba_one),
      'K_a_to_b':str(K_ab_one),
  },
  'symmetric_two_cycle_cases':len(rows),
  'symmetric_max_cos_theta':str(cstar),
  'symmetric_max_edge_weight':str(fstar),
  'expected_max_2_over_3sqrt3':str(expected),
  'max_derivative_residual':str(fprime),
  'sign_taxonomy_cases':sign_cases,
  'interpretation':(
      'The exact angular Biot-Savart stretching kernel is genuinely directed: one-way edges exist and opposite transactions are not antisymmetric. '
      'A symmetric 3D configuration supports a positive two-cycle in which both vorticity elements stretch each other, with maximal dimensionless angular weight 2/(3*sqrt(3)) at cos(theta)=1/sqrt(3). '
      'Same-sign two-cycles require opposite longitudinal projections along the separation axis and a nonzero triple product. '
      'Therefore vortex-stretching transactions cannot be booked as pairwise transfer of a conserved scalar; the next obstruction must concern the lifetime/geometry of productive cycles.'
  ),
  'symmetric_family':rows,
},indent=2,allow_nan=False))
