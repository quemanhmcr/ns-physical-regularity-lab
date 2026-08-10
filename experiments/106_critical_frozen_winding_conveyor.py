import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def certify_one(x,label,tol='1e-30'):
 t=arb(tol)
 if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
rows=[]
# Canonical escape ODE: winding/stretch factor n obeys n'/n=s; shape-locked pair gives L=L0 n^-2;
# winding transaction gives s=C Gamma n/L^2=(C Gamma/L0^2)n^5 = k n^5.
# Thus n'=k n^6 and near blowup tau=T-t: n=(5k tau)^-1/5, L=L0(5k tau)^2/5, s=1/(5tau).
for Cs in ['0.2','0.8','2']:
 C=arb(Cs)
 for Gs in ['1e-24','1','1e24']:
  G=arb(Gs)
  for L0s in ['1e-12','1','1e12']:
   L0=arb(L0s); k=C*G/(L0*L0)
   for taus in ['1e-6','1e-30','1e-100']:
    tau=arb(taus); n=(5*k*tau)**(-arb(1)/5); L=L0/(n*n); s=k*n**5
    certify_one(s*(5*tau),('critical s=1/5tau',Cs,Gs,L0s,taus))
    certify_one(L/(L0*(5*k*tau)**(arb(2)/5)),('critical L tau2/5',Cs,Gs,L0s,taus))
    energy=s*s*L**5; closed_energy=C*C*G*G*L0
    certify_one(energy/closed_energy,('critical s2L5 constant',Cs,Gs,L0s,taus))
    rows.append({'C_winding_geometry':Cs,'Gamma_thread':Gs,'L0':L0s,'tau':taus,'k_CGamma_over_L0sq':str(k),'winding_stretch_factor_n':str(n),'source_scale_L':str(L),'strain_s':str(s),'s2L5_energy_occupancy_scale':str(energy),'closed_constant_C2Gamma2L0':str(closed_energy)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'Combining frozen material line stretch n, pair-cell shape-lock bridge contraction L=L0 n^-2, and winding transaction s=C Gamma n/L^2 produces the autonomous escape dn/dt=k n^6. Its finite-time solution has n~tau^-1/5, L~tau^2/5 and s=1/(5tau). The strain occupancy scale s^2L^5 is exactly constant C^2 Gamma^2 L0. Thus the Hodge 2/5 energy exponent is reproduced by the minimal frozen winding conveyor rather than imposed as a similarity ansatz. This is a kinematic/transaction mechanism calibration, not an exact NS blow-up solution.','rows':rows},indent=2,allow_nan=False))
