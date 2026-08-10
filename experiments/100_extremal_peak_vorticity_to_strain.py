import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def root(x,p): return (x.log()/p).exp()
rows=[]
# E=diag(2,-1,-1): max_n |n cross E n|=3/2, hence max angular |omega_prod|=(5/2)q.
for p_int in [1,2,4,8,32,128]:
 p=arb(p_int); Cp=5*p*(p+9)/(14*(p+2)*(p+7)); xm=root(arb(2)/(p+2),p)
 qmax_over_q0=xm*xm*p/(p+2)
 qmax_over_s=qmax_over_q0/Cp
 omega_peak_over_s=(arb(5)/2)*qmax_over_s
 if not (omega_peak_over_s>arb(2) and omega_peak_over_s<arb(7)): raise AssertionError(('peak productive vorticity not O(strain)',p_int,omega_peak_over_s))
 rows.append({'p':p_int,'Hodge_strain_coefficient_Cp':str(Cp),'radial_peak_xm':str(xm),'qmax_over_q0':str(qmax_over_q0),'qmax_over_source_strain':str(qmax_over_s),'peak_minimum_productive_vorticity_over_source_strain':str(omega_peak_over_s)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'For the localized sharp extremal carrier q=q0 x^2(1-x^p), the productive shell amplitude peaks at x_m=[2/(p+2)]^(1/p). With E=diag(2,-1,-1), max angular |n cross E n|=3/2, so the minimum productive vorticity peak is (5/2)q_max. Its ratio to the Hodge source strain is an order-one number (between 2 and 7 in the tested family, tending to 7 as p grows). Thus high-Re source recruitment cannot be supplied by weak background fluid in this extremal geometry: the material entering the productive source must already carry vorticity on the current strain scale.','rows':rows},indent=2,allow_nan=False))
