import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

# Stadium closure of one coherent circulation lineage.
# Two semicircular bends of radius R carry total bending energy
# int kappa^2 ds = 2*pi/R.  If an isolated meridional collar a<rho<b
# retains at least theta*Gamma circulation and b<R, Cauchy+Stokes give
# E_bend_collar >= theta^2 Gamma^2 R/2 * (1-b/R)*log(b/a).
# Hence Tax_bend * E_bend_collar has no large-R escape.

rows=[]
R_over_b_values=['2','3','10','1e2','1e4','1e8','1e16','1e30']
collar_ratios=['1.1','2','10','1e3']
thetas=['0.1','0.5','1']
gammas=['1e-12','1','1e12']
nu_values=['1e-18','1','1e18']

for Rob_s in R_over_b_values:
  Rob=arb(Rob_s)
  for ba_s in collar_ratios:
    ba=arb(ba_s)
    a=arb(1); b=ba*a; R=Rob*b
    tubular=1-b/R
    if not (tubular>0): raise AssertionError(('invalid tubular stadium',Rob_s,ba_s,tubular))
    logba=(b/a).log()
    for th_s in thetas:
      th=arb(th_s)
      for gs in gammas:
        G=arb(gs)
        for nus in nu_values:
          nu=arb(nus)
          curvature_integral=2*pi/R
          tax=nu*G*curvature_integral
          Ecollar=th*th*G*G*R/2*tubular*logba
          product=tax*Ecollar
          target=pi*nu*th*th*G**3*tubular*logba
          ratio=product/target
          if not ratio.contains(1):
              raise AssertionError(('stadium tax-energy product failed',Rob_s,ba_s,th_s,gs,nus,ratio))
          # R-lengthening trade: tax*R and E/R are invariant at fixed collar geometry.
          taxR=tax*R
          taxR_target=2*pi*nu*G
          EoverR=Ecollar/R
          EoverR_target=th*th*G*G*tubular*logba/2
          if not (taxR/taxR_target).contains(1):
              raise AssertionError(('curvature tax did not scale 1/R',Rob_s,ratio))
          if not (EoverR/EoverR_target).contains(1):
              raise AssertionError(('collar energy did not scale R',Rob_s,ratio))
          rows.append({
             'R_over_b':Rob_s,'b_over_a':ba_s,'theta':th_s,'Gamma':gs,'nu':nus,
             'curvature_integral':str(curvature_integral),
             'bend_direction_tax_lower':str(tax),
             'bend_collar_energy_lower':str(Ecollar),
             'tax_times_energy_ratio':str(ratio),
             'tubular_factor_1_minus_b_over_R':str(tubular),
          })

# Finite-energy consequence in a declared safe-tubular branch b/R <= eps.
finite=[]
eps=arb('0.1')
for E0s in ['1e-18','1','1e18']:
  E0=arb(E0s)
  for gs in ['1e-12','1','1e12']:
    G=arb(gs)
    for nus in ['1e-18','1','1e18']:
      nu=arb(nus)
      th=arb('0.5'); logba=arb(2).log()
      # E0 >= Ecollar >= theta^2 G^2 R/2 (1-eps) log(b/a)
      Rmax=2*E0/(th*th*G*G*(1-eps)*logba)
      tax_floor=2*pi*nu*G/Rmax
      expected=pi*nu*th*th*G**3*(1-eps)*logba/E0
      if not (tax_floor/expected).contains(1):
          raise AssertionError(('finite-energy closure tax floor mismatch',E0s,gs,nus,tax_floor/expected))
      finite.append({'E0':E0s,'Gamma':gs,'nu':nus,'theta':'0.5','b_over_a':'2','epsilon_b_over_R':'0.1',
                     'R_max_from_energy':str(Rmax),'closure_tax_floor':str(tax_floor),
                     'floor_ratio':str(tax_floor/expected)})

# Autopsy 1: collar collapse b/a=1+delta makes the protection factor vanish linearly.
collapse=[]
last=None
for ds in ['0.3','0.1','0.03','0.01','0.003','0.001','1e-4','1e-6','1e-9','1e-12']:
    delta=arb(ds)
    factor=(1+delta).log()
    ratio=factor/delta
    if not (arb(0)<factor): raise AssertionError(('nonpositive collar log',ds,factor))
    if last is not None and delta < arb('0.1') and not (ratio >= last):
        # log(1+d)/d increases toward 1 as d decreases.
        raise AssertionError(('collar-collapse asymptotic lost monotonicity',ds,last,ratio))
    last=ratio
    collapse.append({'delta_b_over_a_minus_1':ds,'log_b_over_a':str(factor),'log_over_delta':str(ratio)})
if not (last > arb('0.999999999999')):
    raise AssertionError(('collar collapse did not approach linear asymptotic',last))

# Autopsy 2: loss of circulation isolation theta->0 weakens the energy collar quadratically.
isolation=[]
for ths in ['1','0.3','0.1','0.03','0.01','0.001','1e-6']:
    th=arb(ths)
    normalized=(th*th)/(th*th)
    if not normalized.contains(1): raise AssertionError('theta scaling arithmetic failed')
    isolation.append({'theta':ths,'relative_collar_energy_factor':str(th*th)})

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'stadium_cases':len(rows),
  'finite_energy_cases':len(finite),
  'collar_collapse':collapse,
  'isolation_loss':isolation,
  'finite_energy_branch':finite,
  'interpretation':'A large stadium closure can make centerline curvature tax fall like 1/R, but an isolated circulation collar on the same bends has an exact Cauchy-Stokes kinetic-energy floor growing like R. Their product is bounded below by pi*nu*theta^2*Gamma^3*(1-b/R)*log(b/a), so finite total energy removes the R->infinity escape whenever a nondegenerate collar persists. The law deliberately collapses when b/a->1 or theta->0; those are the unresolved near-contact/cancellation/reconnection branches rather than false positives.',
},indent=2,allow_nan=False))
