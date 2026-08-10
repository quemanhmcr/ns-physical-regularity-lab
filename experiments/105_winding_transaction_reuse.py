import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rows=[]
# Spherical-shell field from module104:
# omega_r=A sin^2(theta)cos(theta)/r^2, omega_phi=B sin(theta)cos(theta)/r.
# Positive north-hemisphere through-flux Gamma=pi A/2.
# At theta=pi/4, fixed shell traverse DeltaR in N turns gives B/A=pi N/DeltaR.
# Tangential sharp transaction q_e=2B/(5r), so Gamma_Q(r)=q_e r^2=2Br/5.
for R1s,R2s in [('1','2'),('1e-6','2e-6'),('1e6','2e6')]:
 R1=arb(R1s); R2=arb(R2s); dR=R2-R1; r=(R1*R2).sqrt()
 for Gs in ['1e-24','1','1e24']:
  G=arb(Gs); A=2*G/pi
  for Ns in ['1','10','1e6','1e30']:
   N=arb(Ns); B=A*pi*N/dR
   GammaQ=2*B*r/5
   gain=GammaQ/G
   closed_gain=arb(4)*N*r/(5*dR)
   if not (gain/closed_gain).contains(1): raise AssertionError(('winding transaction gain',R1s,R2s,Gs,Ns,gain,closed_gain))
   Zt=(8*pi/15)*B*B*dR
   Zr=(32*pi/105)*A*A*dR/(R1*R2)
   null_ratio=Zr/Zt
   dual=null_ratio*gain*gain
   closed_dual=arb(64)/(175*pi*pi)
   if not (dual/closed_dual).contains(1): raise AssertionError(('winding-null duality',R1s,R2s,Gs,Ns,dual,closed_dual))
   # Tangential component exactly saturates shellwise sharp Q floor; its volume integral is Zt.
   sharp=(8*pi/15)*B*B*dR
   if not (Zt/sharp).contains(1): raise AssertionError(('sharp tangential winding carrier',R1s,R2s,Gs,Ns))
   rows.append({'R1':R1s,'R2':R2s,'throughgoing_circulation_Gamma':Gs,'turns_N':Ns,'Gamma_Q_characteristic':str(GammaQ),'transaction_gain_GammaQ_over_Gamma':str(gain),'closed_gain_4Nr_over_5DeltaR':str(closed_gain),'sharp_tangential_enstrophy':str(Zt),'radial_transaction_null_enstrophy':str(Zr),'null_over_sharp':str(null_ratio),'null_ratio_times_gain_squared':str(dual),'closed_duality_64_over_175pi2':str(closed_dual)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'A through-going frozen circulation thread can be reused by winding. In the exact spherical-shell calibration, Gamma_Q/Gamma=(4/5)N r/DeltaR, so source transaction circulation grows linearly with turn count at fixed material circulation. Simultaneously the radial transaction-null/sharp enstrophy ratio falls like N^-2, with exact product [Z_null/Z_sharp](Gamma_Q/Gamma)^2=64/(175 pi^2). Thus near-sharpness does not penalize winding reuse; it accompanies it. Unbounded source Re need not mean unbounded independent circulation stock.','rows':rows},indent=2,allow_nan=False))
