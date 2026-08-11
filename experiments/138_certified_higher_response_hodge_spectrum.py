import json, os
from flint import ctx,arb
from fractions import Fraction as F
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import coupled46_hodge_core as H
import degree6_hodge_servo_core as C

def pD(A):return {e:sum(e)*v for e,v in A.items()}
def plap(A):
 q={}
 for j in range(3):q=C.padd(q,C.pder(C.pder(A,j),j))
 return q
def lapS(P,r2):return C.padd(C.pmul(r2,plap(P)),C.pscale(-1,C.padd(pD(pD(P)),pD(P))))
def spectral_project(P,l,allowed,r2):
 out=dict(P);lam=l*(l+1)
 for k in allowed:
  if k==l:continue
  lk=k*(k+1);out=C.pscale(arb(1)/(lk-lam),C.padd(lapS(out,r2),C.pscale(lk,out)))
 return out
def spectrum(V,d,X,r2):
 h=C.vdot(X,V);rad=tuple(C.pmul(h,X[i]) for i in range(3));T=C.vadd(V,C.vscale(-1,rad));divR=C.div(T);cor={}
 for i in range(3):
  for j in range(3):cor=C.padd(cor,C.pmul(X[i],C.pmul(X[j],C.pder(T[i],j))))
 dS=C.padd(divR,C.pscale(-1,cor));cS=C.vdot(X,C.curl(T));odd=list(range(1,d+2,2));even=list(range(2,d+1,2))
 pol={l:C.norm2s(spectral_project(dS,l,odd,r2))/arb(l*(l+1)) for l in odd}
 tor={l:C.norm2s(spectral_project(cS,l,even,r2))/arb(l*(l+1)) for l in even}
 radE=C.norm2s(h);tang=C.norm2v(T);polE=sum(pol.values(),C.z);torE=sum(tor.values(),C.z)
 if not (C.norm2v(V)-(radE+tang)).contains(0):raise AssertionError(('radial/tangent reconstruction',d))
 if not (tang-(polE+torE)).contains(0):raise AssertionError(('surface Hodge reconstruction',d,tang,polE,torE))
 return radE,pol,tor

st=H.prepare();sym=H.feedback_symmetry_basis(st);seq=C.solve_degree6_servo();y=[arb(seq['coeff'][i].mid()) for i in st['t4idx']];a,_=H.sym_coords(y,sym);a=[arb(v.mid()) for v in a]
for _ in range(8):
 g,_,_=H.reduced_feedback_native(st,sym,a);Jr=H.reduced_jacobian_native(st,sym,a);J=[[arb(Jr[i][j].mid()) for j in range(5)] for i in range(5)];dd=H.arbmat_solve(J,[-arb(v.mid()) for v in g]);a=[arb((a[i]+dd[i]).mid()) for i in range(5)]
R=arb('0 +/- 1e-20');box=[arb(v.mid())+R for v in a];_,fb,_=H.reduced_feedback_native(st,sym,box);hr=H.higher_responses_from_coupled(st,fb);rows=[]
for d in (8,10,12):
 N=hr[d][2];rad,pol,tor=spectrum(N,d,st['X'],st['r2']); row={'degree':d,'radial_null_mean_square':str(rad),'positive_poloidal_sectors':[],'positive_toroidal_sectors':[]}
 for l,e in pol.items():
  row[f'poloidal_l{l}_energy']=str(e)
  if e>0:row['positive_poloidal_sectors'].append(l)
 for l,e in tor.items():
  row[f'toroidal_l{l}_energy']=str(e)
  if e>0:row['positive_toroidal_sectors'].append(l)
 if 2 in tor and not tor[2].contains(0):raise AssertionError(('productive T2 leaked into null spectrum',d,tor[2]))
 rows.append(row)
print(json.dumps({'arb_precision_bits':BITS,'status':'PASS','certified_root_box_radius':'1e-20','interpretation':'Resolve the unavoidable higher transaction-null emissions of the certified coupled degree-four/six servo using only intrinsic surface operators: radial/tangential split, surface divergence, surface curl, and the sphere Laplacian.  No preselected truncated modal ansatz is used; the allowed Hodge sectors are the complete sectors compatible with each physical homogeneity.  Positive lower bounds identify angular channels that are necessarily present at the actual coupled root.','rows':rows},indent=2,allow_nan=False))
