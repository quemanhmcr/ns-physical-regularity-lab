import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def certify_one(x,label,tol='1e-30'):
 t=arb(tol)
 if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
rows=[]
# Exact calibration: omega=B(t)e_z on a disk of radius R(t), u=0.
# Choose j=(kappa/2)(-y,x,0), so curl j=kappa e_z and Bdot=-kappa satisfies partial_t Omega=-d j.
# Boundary velocity V=Rdot e_r. Then recruitment integral=2pi B R Rdot and int_boundary j.dl=pi kappa R^2.
for Bs in ['1e-24','1','1e24']:
 B=arb(Bs)
 for Rs in ['1e-12','1','1e12']:
  R=arb(Rs)
  for vs in ['-1e6','-0.1','0','0.1','1e6']:
   v=arb(vs)
   for ks in ['1e-24','1','1e24']:
    k=arb(ks)
    flux=pi*R*R*B
    fluxdot=2*pi*B*R*v-pi*k*R*R
    recruit=2*pi*B*R*v
    viscous=pi*k*R*R
    if viscous==0:
      if not (fluxdot-recruit).contains(0): raise AssertionError(('moving-surface ledger k=0',Bs,Rs,vs))
    else:
      certify_one((recruit-viscous)/fluxdot,('moving-source flux ledger',Bs,Rs,vs,ks)) if not fluxdot.contains(0) else None
      if not (fluxdot-(recruit-viscous)).contains(0): raise AssertionError(('exact flux ledger',Bs,Rs,vs,ks,fluxdot,recruit,viscous))
    rows.append({'B':Bs,'R':Rs,'Rdot':vs,'kappa':ks,'source_flux_piR2B':str(flux),'direct_flux_derivative':str(fluxdot),'material_recruitment_boundary_integral':str(recruit),'Kelvin_viscous_boundary_integral':str(viscous),'ledger_recruitment_minus_viscous':str(recruit-viscous)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For a non-material moving source surface Sigma(t), the exact vorticity-flux transport law is d/dt int_Sigma Omega = int_boundary [i_{V-u}Omega - j]. The calibration combines a changing uniform vorticity flux, a moving disk boundary, and a nonzero Kelvin current. It separates material recruitment caused by source-boundary motion from viscous ancestry replacement. A shrinking Eulerian source can therefore change its ancestry inventory even in the high-Re approximately frozen branch by cutting through material flux.','rows':rows},indent=2,allow_nan=False))
