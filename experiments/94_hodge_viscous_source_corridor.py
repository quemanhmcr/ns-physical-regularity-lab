import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def root(x,n): return (x.log()/n).exp()
def certify_one(x,label,tol='1e-30'):
 t=arb(tol)
 if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))
rows=[]
for Es in ['1e-24','1','1e24']:
 E=arb(Es)
 for nus in ['1e-12','1','1e12']:
  nu=arb(nus)
  sc=(arb(28)**5*pi*pi/(arb(900)))*nu**5/(E*E)
  for mults in ['1e-12','1','1e12']:
   mult=arb(mults); s=sc*mult
   LE=root(arb(30)*E/(pi*s*s),5)
   Lv=(arb(28)*nu/s).sqrt()
   ratio=LE/Lv
   tenth=ratio**10
   certify_one(tenth/(s/sc),('corridor tenth-power identity',Es,nus,mults))
   ReE=s*LE*LE/nu
   fifth=root(s/sc,5)
   certify_one(ReE/(arb(28)*fifth),('energy-horizon Re identity',Es,nus,mults))
   rows.append({'E0':Es,'nu':nus,'strain_over_critical':mults,'critical_strain_sc':str(sc),'strain_s':str(s),'Hodge_energy_horizon_LE':str(LE),'unit_exposure_viscous_length_Lnu':str(Lv),'corridor_ratio_LE_over_Lnu':str(ratio),'ratio_tenth':str(tenth),'Re_at_energy_horizon':str(ReE),'corridor_nonempty':bool(LE>=Lv)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'Combining the half-strain Hodge energy horizon LE=[30E0/(pi s^2)]^(1/5) with the localized spectral-gap unit-exposure length Lnu=sqrt(28nu/s) gives an intrinsic source corridor. Its exact identities are (LE/Lnu)^10=s/sc and Re(LE)=28(s/sc)^(1/5), where sc=(28^5 pi^2/900)nu^5/E0^2. Above sc a high-Re source can live between the viscous length and the energy-forced Hodge horizon; below sc the unit-exposure viscous length is already outside the energy horizon.','rows':rows},indent=2,allow_nan=False))
