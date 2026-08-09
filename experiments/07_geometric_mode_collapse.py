import json, os, math
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def dot(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def add(a,b): return (a[0]+b[0],a[1]+b[1],a[2]+b[2])
def sub(a,b): return (a[0]-b[0],a[1]-b[1],a[2]-b[2])
def mul(c,a): return (c*a[0],c*a[1],c*a[2])

def loop(th,eps,m):
 mm=arb(m); c=th.cos(); s=th.sin()
 return ((c,s,eps*(mm*th).sin()),(-s,c,eps*mm*(mm*th).cos()))

def sigma(theta,eps,m,core,N):
 x,xt=loop(theta,eps,m); h=2*pi/N; su=(arb(0),arb(0),arb(0))
 for j in range(N):
  ph=h*j; xp,xpp=loop(ph,eps,m); r=sub(x,xp)
  d2=dot(r,r)+core*core; d32=d2*d2.sqrt(); d52=d32*d2
  term=add(mul(1/d32,cross(xpp,xt)),mul(-3*dot(r,xt)/d52,cross(xpp,r)))
  su=add(su,mul(h,term))
 uth=mul(1/(4*pi),su)
 return dot(xt,uth)/dot(xt,xt)

def eff(eps,m,core,N,probes=17):
 ss=arb(0)
 for k in range(probes):
  th=2*pi*(arb(k)+arb('0.193'))/probes
  s=sigma(th,eps,m,core,N); ss+=s*s
 rms=(ss/probes).sqrt()
 return 4*pi*rms  # R=Gamma=1

eps=arb('1e-5')
modes=[2,3,4,6,8]
cores=['0.06','0.12','0.24']
rows=[]
for c_s in cores:
 c=arb(c_s)
 for m in modes:
  e256=eff(eps,m,c,256)
  e512=eff(eps,m,c,512)
  coeff=e512/eps
  geom=coeff/(m*m-1)
  rel=abs(float(e512.mid())-float(e256.mid()))/max(abs(float(e512.mid())),1e-300)
  if rel>0.12:
   raise AssertionError(f'quadrature unstable m={m} c={c_s}: {rel}')
  rows.append({'mode':m,'core_over_R':c_s,'q=m*core':str(arb(m)*c),'eff':str(e512),'coeff_eff_over_eps':str(coeff),'normalized_by_m2minus1':str(geom),'refinement_rel_mid':rel,'norm_mid':float(geom.mid())})

# Test the stronger collapse claim at matched q. We do not force PASS: report its spread and classify.
matched={
 '0.24':[(2,'0.12'),(4,'0.06')],
 '0.48':[(2,'0.24'),(4,'0.12'),(8,'0.06')],
 '0.72':[(3,'0.24'),(6,'0.12')],
}
lookup={(r['mode'],r['core_over_R']):r['norm_mid'] for r in rows}
collapse=[]
max_spread=0.0
for q,pairs in matched.items():
 vals=[lookup[p] for p in pairs]
 mean=sum(vals)/len(vals); spread=(max(vals)-min(vals))/max(abs(mean),1e-300)
 max_spread=max(max_spread,spread)
 collapse.append({'q':q,'pairs':pairs,'normalized_values':vals,'relative_spread':spread})

classification='SUPPORTED' if max_spread<0.08 else ('APPROXIMATE' if max_spread<0.25 else 'REJECTED')
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','eps':'1e-5','cases':len(rows),
 'mode_factor_hypothesis':'eff ~ eps*(m^2-1)*F(m*core/R)',
 'matched_q_collapse_classification':classification,'max_matched_q_relative_spread':max_spread,
 'interpretation':'The exact planar zero mode requires a geometric factor vanishing at m=1. This experiment tests whether the leading 3D self-stretching coefficient factors into (m^2-1) times a core-scale transfer function; the classification is diagnostic, not assumed.',
 'rows':rows,'matched_q_groups':collapse
},indent=2))
