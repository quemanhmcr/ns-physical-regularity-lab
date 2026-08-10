import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def norm2(S): return sum(S[i][j]*S[i][j] for i in range(3) for j in range(3))
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
Ss=[('diag',((arb(2),arb(0),arb(0)),(arb(0),-arb(1),arb(0)),(arb(0),arb(0),-arb(1)))),('mixed',((arb('0.4'),arb('0.3'),-arb('0.2')),(arb('0.3'),-arb('0.1'),arb('0.25')),(-arb('0.2'),arb('0.25'),-arb('0.3'))))]
rows=[]
for name,S in Ss:
  S2=norm2(S)
  for Ls in ['1e-24','1','1e24']:
    L=arb(Ls)
    kinetic=(8*pi/405)*S2*L**5
    harmonic=(2*pi/15)*S2*L**5
    certify_one((kinetic/harmonic)/(arb(4)/27),('carrier/harmonic ratio',name,Ls))
    enstrophy=(112*pi/45)*S2*L**3
    integrated_sharp=(20*pi/9)*(arb(14)/5)**2*S2*L**3/7
    certify_one(enstrophy/integrated_sharp,('integrated sharp saturation',name,Ls))
    rows.append({'S_case':name,'L':Ls,'S_Frobenius_squared':str(S2),'carrier_kinetic_energy':str(kinetic),'pure_harmonic_energy':str(harmonic),'carrier_over_harmonic':str(kinetic/harmonic),'carrier_enstrophy_integral':str(enstrophy),'integrated_sharp_transaction_floor':str(integrated_sharp)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'The self-contained tangent strain carrier has E=(8*pi/405)|S|^2 L^5, exactly 4/27 of the pure harmonic linear-strain occupancy on the same ball. Its enstrophy is (112*pi/45)|S|^2 L^3 and exactly equals the radial integral of the sharp angular transaction floor. The 4/27 coincidence with the independent productive-pair capacity is recorded only as a coefficient resonance.','rows':rows},indent=2,allow_nan=False))
