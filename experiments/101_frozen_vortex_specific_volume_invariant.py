import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def mv(A,v): return tuple(sum(A[i][j]*v[j] for j in range(3)) for i in range(3))
def norm(v): return dot(v,v).sqrt()
def certify_one(x,label,tol='1e-30'):
 t=arb(tol)
 if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
rows=[]
omegas=[(arb(1),arb(2),arb(3)),(arb('-0.3'),arb('0.7'),arb('1.1'))]
for ls in ['1e-24','1','1e24']:
 lam=arb(ls)
 for ks in ['-1e24','0','1','1e24']:
  k=arb(ks); F=((lam,k,arb(0)),(arb(0),1/lam,arb(0)),(arb(0),arb(0),arb(1)))
  for w0 in omegas:
   m0=norm(w0); xi0=tuple(x/m0 for x in w0); ds0=arb('0.7')
   wt=mv(F,w0); dxt=mv(F,xi0); dst=norm(dxt)*ds0; mt=norm(wt)
   certify_one((dst/mt)/(ds0/m0),('frozen ds-over-omega invariant',ls,ks,w0))
   # Flux-tube infinitesimal identity dV=dGamma ds/|omega|.
   dGamma=arb('0.13'); dV0=dGamma*ds0/m0; dVt=dGamma*dst/mt
   certify_one(dVt/dV0,('specific volume per circulation invariant',ls,ks,w0))
   rows.append({'lambda':ls,'shear_k':ks,'omega0':[str(x) for x in w0],'ds0':str(ds0),'omega0_mag':str(m0),'dst':str(dst),'omegat_mag':str(mt),'ds_over_omega_ratio':str((dst/mt)/(ds0/m0)),'specific_flux_volume_ratio':str(dVt/dV0)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For incompressible frozen-vorticity Cauchy transport omega=F omega0, a material line element initially tangent to omega obeys ds(t)/|omega(t)|=ds0/|omega0| pointwise, even under extreme shear and anisotropic stretch. For a vortex-flux element dGamma this is exactly dV=dGamma ds/|omega|, so mu=int ds/|omega| is the material specific volume per circulation flux of a closed lineage.','rows':rows},indent=2,allow_nan=False))
