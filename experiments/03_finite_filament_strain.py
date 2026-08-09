import json, os
from flint import arb, ctx

BITS = int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec = BITS
pi = arb.pi()

alpha_values = ['1e-12','1e-8','1e-4','0.01','0.1','0.5','1','2','10','1e4','1e8','1e12']
d_values = ['1e-12','1e-3','1','1e6','1e12']
gamma_values = ['1e-12','1','1e12']
rows=[]
for a_s in alpha_values:
    alpha = arb(a_s)
    formula = alpha*(alpha*alpha+2)/(1+alpha*alpha)**arb('1.5')
    fprime = (2-alpha*alpha)/(1+alpha*alpha)**arb('2.5')
    for d_s in d_values:
        d=arb(d_s)
        L=2*alpha*d
        for g_s in gamma_values:
            gamma=arb(g_s)
            root=(d*d+(L/2)*(L/2)).sqrt()
            uy=gamma*L/(4*pi*d*root)
            du_dd = -gamma*L/(4*pi)*(1/(d*d*root)+1/(root*root*root))
            sxy = -du_dd/2
            efficiency = 4*pi*d*d*sxy/gamma
            ratio = efficiency/formula
            if not ratio.contains(1):
                raise AssertionError(f'efficiency mismatch alpha={a_s},d={d_s},G={g_s}: {ratio}')
            rows.append({'alpha':a_s,'d':d_s,'Gamma':g_s,'efficiency':str(efficiency),'formula':str(formula),'ratio':str(ratio),'fprime':str(fprime),'u_y':str(uy)})

sqrt2=arb(2).sqrt()
fmax=sqrt2*(sqrt2*sqrt2+2)/(1+sqrt2*sqrt2)**arb('1.5')
expected_fmax=4*arb(2).sqrt()/(3*arb(3).sqrt())
if not (fmax/expected_fmax).contains(1):
    raise AssertionError('maximum efficiency identity failed')

e1=arb(3)/(2*arb(2).sqrt())
formula1=arb(1)*(arb(1)+2)/(arb(2)**arb('1.5'))
if not (formula1/e1).contains(1):
    raise AssertionError('alpha=1 exact value mismatch')

print(json.dumps({
  'arb_precision_bits':BITS,'cases':len(rows),'status':'PASS','fmax':str(fmax),'efficiency_alpha_1':str(e1),
  'interpretation':'A finite open vortex segment with L comparable to d already supplies O(Gamma/d^2) strain. A naive scale-independent donor-length lower bound is therefore false; the next physical constraint must use divergence-free vortex-line closure/global geometry.',
  'rows':rows
}, indent=2))
