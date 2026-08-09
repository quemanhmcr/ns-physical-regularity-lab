import json, os, math
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def cross(a,b):
 return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def dot(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def add(a,b): return (a[0]+b[0],a[1]+b[1],a[2]+b[2])
def sub(a,b): return (a[0]-b[0],a[1]-b[1],a[2]-b[2])
def mul(c,a): return (c*a[0],c*a[1],c*a[2])

def loop(th,R,eps,m):
 c,s=th.cos(),th.sin(); mm=arb(m)
 return ((R*c,R*s,eps*R*(mm*th).sin()),
         (-R*s,R*c,eps*R*mm*(mm*th).cos()))

def sigma_at(theta,R,G,eps,m,core,N):
 x,xt=loop(theta,R,eps,m)
 h=2*pi/N
 su=(arb(0),arb(0),arb(0))
 for j in range(N):
  ph=h*j
  xp,xpp=loop(ph,R,eps,m)
  r=sub(x,xp)
  d2=dot(r,r)+core*core
  d32=d2*d2.sqrt(); d52=d32*d2
  term=add(mul(1/d32,cross(xpp,xt)),mul(-3*dot(r,xt)/d52,cross(xpp,r)))
  su=add(su,mul(h,term))
 uth=mul(G/(4*pi),su)
 return dot(xt,uth)/dot(xt,xt)

def rms_eff(R,G,eps,m,core_frac,N,probes=16):
 core=core_frac*R
 ss=arb(0)
 vals=[]
 for k in range(probes):
  th=2*pi*(arb(k)+arb('0.217'))/probes
  sig=sigma_at(th,R,G,eps,m,core,N)
  vals.append(sig)
  ss += sig*sig
 rms=(ss/probes).sqrt()
 eff=4*pi*R*R*rms/G
 return eff, vals

core_frac=arb('0.12')
eps_list=['1e-4','3e-4','1e-3','3e-3','1e-2','3e-2','0.1','0.3']
rows=[]
# null control: m=1 is a tilted planar loop for any eps.
for e_s in ['0.001','0.1','0.5']:
 e=arb(e_s)
 eff,vals=rms_eff(arb(1),arb(1),e,1,core_frac,128)
 for v in vals:
  if not v.contains(0):
   raise AssertionError(f'm=1 planar control produced stretching: eps={e_s}, sigma={v}')
 rows.append({'mode':1,'eps':e_s,'eff_rms':str(eff),'control':'planar-null'})

# genuinely 3D modes, with N refinement to make sure the observed signal is not quadrature noise.
for m in [2,3]:
 for e_s in eps_list:
  e=arb(e_s)
  eff128,_=rms_eff(arb(1),arb(1),e,m,core_frac,128)
  eff256,_=rms_eff(arb(1),arb(1),e,m,core_frac,256)
  mid128=float(eff128.mid()); mid256=float(eff256.mid())
  rel=abs(mid256-mid128)/max(abs(mid256),1e-300)
  if e_s in ['0.01','0.1'] and rel>0.08:
   raise AssertionError(f'quadrature refinement unstable: m={m}, eps={e_s}, rel={rel}')
  ratio=eff256/e
  rows.append({'mode':m,'eps':e_s,'eff_rms_N128':str(eff128),'eff_rms_N256':str(eff256),'eff_over_eps_N256':str(ratio),'refinement_rel_mid':rel})

# NS scaling check: same dimensionless geometry, R and Gamma changed together only through the expected Gamma/R^2 factor.
scale_rows=[]
for R_s,G_s in [('1','1'),('1e-3','1e6'),('1e3','1e-6')]:
 R,G=arb(R_s),arb(G_s)
 eff,_=rms_eff(R,G,arb('0.03'),2,core_frac,192)
 scale_rows.append({'R':R_s,'Gamma':G_s,'dimensionless_eff':str(eff),'mid':float(eff.mid())})
base=scale_rows[0]['mid']
for row in scale_rows[1:]:
 if abs(row['mid']-base)/max(abs(base),1e-300)>0.02:
  raise AssertionError('dimensionless scaling check failed')

print(json.dumps({
 'arb_precision_bits':BITS,
 'status':'PASS',
 'core_over_R':'0.12',
 'interpretation':'Planarity is an exact geometric gate: m=1 remains non-stretching. Genuine 3D nonplanarity (m=2,3) turns on a dimensionless self-stretching channel; the small-amplitude data expose its leading dependence on nonplanarity without replacing Biot-Savart geometry by absolute-value bounds.',
 'rows':rows,
 'scale_checks':scale_rows
},indent=2))
