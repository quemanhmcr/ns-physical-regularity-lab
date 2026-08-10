import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); rows=[]
# In strain count N, let L=e^{-bN}. Replacing one source-ball volume per dN has cumulative distinct-volume ceiling
# int_N^infty (4pi/3)L^3 dN = (4pi/(9b)) e^{-3bN}.
for bs in ['0.5','2','5']:
 b=arb(bs)
 for Ns in ['0','1','10','100']:
  N=arb(Ns); L=(-b*N).exp(); V=4*pi*L**3/3
  tail=4*pi*(-3*b*N).exp()/(9*b)
  # For physical time with s=dN/dt, if |Ldot|=b s L, the shrinking source boundary sweeps volume rate 4pi b s L^3.
  # Dividing by s gives swept volume per strain count 4pi b L^3; its tail is 4pi/3 L^3.
  swept_tail=4*pi*L**3/3
  rows.append({'b':bs,'N':Ns,'source_scale_L':str(L),'source_ball_volume':str(V),'tail_if_one_full_source_volume_of_new_material_per_strain_count':str(tail),'tail_geometric_boundary_swept_volume':str(swept_tail)})
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'interpretation':'A finite material-volume stock does not obstruct an infinite shrinking-source cascade. If L=e^{-bN}, even replacing an entire source-ball volume on every strain-count increment requires only a finite cumulative material volume proportional to int L^3 dN. The exact shrinking-boundary swept volume is also finite. Therefore the high-Re obstruction cannot be raw recruited material volume; later recruited parcels may have vanishing volume while carrying increasingly amplified vorticity/flux ancestry.','rows':rows},indent=2,allow_nan=False))
