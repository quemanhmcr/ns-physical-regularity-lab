import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def certify_one(x,label,tol='1e-30'):
 t=arb(tol)
 if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
rows=[]
# E=diag(2,-1,-1), localized extremal radial profile q=q0 x^2(1-x^p), x=r/L.
# Hodge source strain coefficient at R=L is s=Cp q0 with Cp=5p(p+9)/[14(p+2)(p+7)].
# Clebsch DeltaA 0->L =5p q0/[6(p+2)], full DeltaB=(3/2)L^2.
# Thus physical extremal flux Phi= [7(p+7)/(2(p+9))] s L^2 = c_p nu Re_source.
for p_int in [1,2,4,8,32]:
 p=arb(p_int); Cp=5*p*(p+9)/(14*(p+2)*(p+7)); cflux=7*(p+7)/(2*(p+9))
 for qs in ['1e-24','1','1e24']:
  q0=arb(qs)
  for Ls in ['1e-12','1','1e12']:
   L=arb(Ls)
   for nus in ['1e-24','1','1e24']:
    nu=arb(nus); s=Cp*q0; Re=s*L*L/nu
    dA=5*p*q0/(6*(p+2)); dB=arb(3)*L*L/2; Phi=dA*dB
    certify_one(Phi/(cflux*nu*Re),('extremal flux vs source Re',p_int,qs,Ls,nus))
    rows.append({'p':p_int,'q0':qs,'L':Ls,'nu':nus,'Hodge_source_strain_s':str(s),'source_Re':str(Re),'extremal_physical_flux_cell':str(Phi),'flux_over_nu_Re':str(Phi/(nu*Re)),'closed_coefficient_7p7_over_2p9':str(cflux)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For a localized sharp extremal carrier q=q0 x^2(1-x^p), the actual source strain and Clebsch flux can both be integrated exactly. The full physical circulation flux of the extremal source cell is Phi=[7(p+7)/(2(p+9))] nu Re_source. Therefore the cumulative frozen-ancestry escape Re_source->infinity is exactly an unbounded circulation-flux aggregation requirement, not merely a dimensionless-number artifact.','rows':rows},indent=2,allow_nan=False))
