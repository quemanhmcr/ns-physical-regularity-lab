import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

# Exact steady Burgers vortex viewed along material circles.
# R_dot=-(a/2)R and q=a R^2/(4 nu), hence q_dot=-a q.
# Material circulation Gamma_m=Gamma_inf(1-exp(-q)) changes even though
# the Eulerian vortex is steady. Kelvin viscosity gives the same derivative.

as_=['1e-30','1','1e30']
nus=['1e-30','1','1e30']
Gs=['1e-18','1','1e18']
q0s=['0.01','0.1','1','3','10','100']
Ns=['0','0.1','1','3','10']  # number of strain times a t
rows=[]
for a_s in as_:
  a=arb(a_s)
  for nu_s in nus:
    nu=arb(nu_s)
    for Gs_ in Gs:
      G=arb(Gs_)
      for q0_s in q0s:
        q0=arb(q0_s)
        for Ns_ in Ns:
          N=arb(Ns_)
          q=q0*(-N).exp()
          emq=(-q).exp()
          Gamma_m=G*(1-emq)
          dGamma=-a*G*q*emq
          # Material radius reconstructed from q.
          R2=4*nu*q/a
          R=R2.sqrt()
          Rdot=-a*R/2
          qdot=-a*q
          # Kelvin current around the moving material circle.
          omega=G*a/(4*pi*nu)*emq
          curlw_theta=omega*a*R/(2*nu)
          kelvin=-nu*(2*pi*R)*curlw_theta
          if not (dGamma/kelvin).contains(1):
              raise AssertionError(('Burgers Kelvin renewal mismatch',a_s,nu_s,Gs_,q0_s,Ns_,dGamma,kelvin))
          if not ((2*Rdot/(R))/(-a)).contains(1):
              raise AssertionError(('material radius contraction mismatch',a_s,nu_s,q0_s,Ns_))
          if not ((qdot/q)/(-a)).contains(1):
              raise AssertionError(('q material clock mismatch',a_s,nu_s,q0_s,Ns_))
          renewal_over_a=(q*emq)/(1-emq)
          renewal_over_transverse=2*renewal_over_a
          rows.append({
            'a':a_s,'nu':nu_s,'Gamma_inf':Gs_,'q0':q0_s,'strain_times_N':Ns_,
            'q_material':str(q),
            'R2':str(R2),
            'material_circulation':str(Gamma_m),
            'dGamma_dt':str(dGamma),
            'Kelvin_viscous_dGamma_dt':str(kelvin),
            'fractional_renewal_rate_over_a':str(renewal_over_a),
            'renewal_over_transverse_compression_rate':str(renewal_over_transverse),
            'ancestry_clock_Chi':str(4*q),
          })

# Rigorous bracket for crossover 2 q e^-q/(1-e^-q)=1, equivalently e^q-1=2q.
ql=arb('1.256')
qh=arb('1.257')
def crossover_fn(q):
    return 2*q*(-q).exp()/(1-(-q).exp())-1
fl=crossover_fn(ql)
fh=crossover_fn(qh)
if not (fl>0 and fh<0):
    raise AssertionError(('Burgers renewal/compression crossover bracket failed',fl,fh))

# The Eulerian profile is steady, but material ancestry changes: certify a concrete case.
q0=arb(10)
q1=q0*(-arb(1)).exp()
f0=1-(-q0).exp(); f1=1-(-q1).exp()
if not (f1<f0):
    raise AssertionError('material circulation did not decay through one strain time')

print(json.dumps({
  'arb_precision_bits':BITS,
  'status':'PASS',
  'cases':len(rows),
  'crossover_q_bracket':['1.256','1.257'],
  'crossover_Chi_bracket':[str(4*ql),str(4*qh)],
  'crossover_function_at_low':str(fl),
  'crossover_function_at_high':str(fh),
  'one_strain_time_q0_10_material_fraction_ratio':str(f1/f0),
  'interpretation':(
      'The exact steady Burgers vortex is an Eulerian structure with continuously renewing material circulation ancestry. '
      'A material circle contracts with q=q0 exp(-a t), while its circulation decays by the Kelvin viscous current. '
      'The fractional renewal rate crosses the transverse material-compression rate near q=1.256.. (Chi about 5.03). '
      'The numerical crossover is model-specific; the promoted mechanism is that steady Eulerian coherence can conceal ongoing material-lineage replacement.'
  ),
  'rows':rows,
},indent=2,allow_nan=False))
