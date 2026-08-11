import json, os
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
z=arb(0);one=arb(1)

def mm(A,B):return tuple(tuple(sum((A[i][k]*B[k][j] for k in range(3)),z) for j in range(3)) for i in range(3))
def mv(A,v):return tuple(sum((A[i][j]*v[j] for j in range(3)),z) for i in range(3))
def tr(A):return sum((A[i][i] for i in range(3)),z)
def det3(A):
 a,b,c=A[0];d,e,f=A[1];g,h,i=A[2]
 return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
def inv3(A):
 a,b,c=A[0];d,e,f=A[1];g,h,i=A[2];D=det3(A)
 if D.contains(0):raise AssertionError(('singular F',D))
 return (((e*i-f*h)/D,(c*h-b*i)/D,(b*f-c*e)/D),((f*g-d*i)/D,(a*i-c*g)/D,(c*d-a*f)/D),((d*h-e*g)/D,(b*g-a*h)/D,(a*e-b*d)/D))
def mn2(A):return sum((A[i][j]*A[i][j] for i in range(3) for j in range(3)),z)
def vn2(v):return sum((x*x for x in v),z)
def msub(A,B):return tuple(tuple(A[i][j]-B[i][j] for j in range(3)) for i in range(3))
def mpow(A,n):
 R=((one,z,z),(z,one,z),(z,z,one))
 for _ in range(n):R=mm(R,A)
 return R

Fs=[
 ('diag',((arb(2),z,z),(z,arb('0.5'),z),(z,z,one))),
 ('shear',((one,arb(2),arb('-0.5')),(z,one,arb(3)),(z,z,one))),
 ('mixed',((arb(3),arb(2),z),(z,arb('0.5'),arb(4)),(z,z,arb(2)/3))),
]
# mixed determinant=1: 3*0.5*(2/3)=1.
Bcases=[
 ('simple',((arb(2),z,z),(z,-one,z),(z,z,-one)),None),
 ('rank2',((one,z,z),(z,-one,z),(z,z,z)),(z,z,one)),
 ('nilpotent_rank2',((z,one,z),(z,z,one),(z,z,z)),(one,z,z)),
]
rows=[]
for fn,F in Fs:
 if not det3(F).contains(1):raise AssertionError(('detF',fn,det3(F)))
 Fi=inv3(F)
 for bn,B0,k0 in Bcases:
  B=mm(mm(F,B0),Fi)
  for k in (1,2,3):
   if not (tr(mpow(B,k))-tr(mpow(B0,k))).contains(0):raise AssertionError(('trace invariant',fn,bn,k))
  if not (det3(B)-det3(B0)).contains(0):raise AssertionError(('det invariant',fn,bn))
  kernel_err='not_applicable'
  if k0 is not None:
   kt=mv(F,k0);res=mv(B,kt)
   if not vn2(res).contains(0):raise AssertionError(('kernel transport',fn,bn,res))
   kernel_err=str(vn2(res))
  nilpotent3='not_applicable'
  if bn=='nilpotent_rank2':
   n3=mn2(mpow(B,3));n2=mn2(mpow(B,2))
   if not n3.contains(0) or not n2.lower()>0:raise AssertionError(('Jordan nilpotent',fn,n2,n3))
   nilpotent3=str(n3)
  rows.append({'F_case':fn,'B0_case':bn,'det_F':str(det3(F)),'det_B0':str(det3(B0)),'det_B_after':str(det3(B)),'tr_B':str(tr(B)),'tr_B2':str(tr(mpow(B,2))),'tr_B3':str(tr(mpow(B,3))),'kernel_transport_error_square':kernel_err,'nilpotent_B3_error_square':nilpotent3})
print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'At an Euler vorticity zero, the Cauchy formula omega(t,X(a,t))=F(a,t)omega0(a) differentiates to B(t)=F B0 F^{-1} because the derivative-of-F term is multiplied by omega0=0.  The artifact certifies the resulting similarity invariants for simple, rank-two, and nontrivial nilpotent zero germs under several volume-preserving deformation gradients. '
  'All spectral traces and determinants are unchanged; for degenerate germs the kernel direction is transported exactly by F, and Jordan nilpotency is preserved.  Thus rank, eigenvalues, determinant/index data and the tangent kernel of a zero manifold are Euler ancestry labels, not quantities that affine or non-affine Euler strain can amplify. '
  'A degenerate zero manifold may still be geometrically folded by the flow, but its local B spectrum must come from its initial zero-germ ancestry unless viscosity mutates it.  This separates spectral amplification from geometric folding and makes the latter the next independent escape mechanism.'),
 'rows':rows
},indent=2,allow_nan=False))
