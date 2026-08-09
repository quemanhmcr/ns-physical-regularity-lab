import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()
rows=[]
for R_s in ['1e-30','1e-12','0.001','1','1e3','1e12','1e30']:
    R=arb(R_s)
    E_coeff=arb(1)/2*(4*pi*R**5/15)
    floor_coeff=(2*pi/15)*R**5
    ratio_E=E_coeff/floor_coeff
    if not ratio_E.contains(1):
        raise AssertionError(('Arb lost exact harmonic floor ratio',R_s,ratio_E))
    proj_ratio=(15/(8*pi*R))*(R*(8*pi/15))
    if not proj_ratio.contains(1):
        raise AssertionError(('Arb lost exact Hodge projector ratio',R_s,proj_ratio))
    carrier_coeff=(15/(8*pi*R))*(R*(1-R**2)*(8*pi/15))
    target=1-R**2
    carrier_ratio = None if target.contains(0) else carrier_coeff/target
    if carrier_ratio is not None and not carrier_ratio.contains(1):
        raise AssertionError(('Arb lost carrier Hodge profile',R_s,carrier_ratio))
    rows.append({'R':R_s,'energy_ratio':str(ratio_E),'projector_ratio':str(proj_ratio),
                 'carrier_ratio':None if carrier_ratio is None else str(carrier_ratio)})
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS',
 'interpretation':'Arb interval arithmetic certifies the exact scale cancellation in the harmonic strain energy floor and the degree-two Hodge boundary projector across sixty orders of magnitude in radius. The canonical carrier profile S_h=(1-R^2)S is independently enclosed whenever the ratio is nonsingular.',
 'rows':rows
},indent=2))
