import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def outer(a,b): return tuple(tuple(a[i]*b[j] for j in range(3)) for i in range(3))
def madd(*Ms): return tuple(tuple(sum(M[i][j] for M in Ms) for j in range(3)) for i in range(3))
def mscale(c,M): return tuple(tuple(c*M[i][j] for j in range(3)) for i in range(3))
def contract(A,B): return sum(A[i][j]*B[i][j] for i in range(3) for j in range(3))
def norm2(M): return contract(M,M)
def Mg(u,v,g): return madd(outer(u,v),outer(v,u),mscale(-g,outer(u,u)),mscale(-g,outer(v,v)))
def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x>1-t and x<1+t): raise AssertionError((label,x))

# Generic productive triad.
a=(arb('0.8'),arb(0),arb('0.6')); n=(arb(0),arb(0),arb(1)); b=(arb('0.3'),arb('-0.4'),-arb(3).sqrt()/2)
coords=[('alpha',a,n,dot(a,n)),('beta',b,n,dot(b,n)),('gamma',a,b,dot(a,b))]
rows=[]
for label,u,v,g in coords:
    M=Mg(u,v,g); M2=norm2(M); closed=2*(1-g*g)**2
    certify_one(M2/closed,('shape tensor norm identity',label))
    for js in ['1e-24','1','1e24']:
      j=arb(js)
      for rs in ['1e-6','1','1e6']:
        r=arb(rs)
        # Minimum-norm harmonic strain producing M:S_h=j is S_h=(j/|M|^2)M.
        Sh=mscale(j/M2,M)
        produced=contract(M,Sh)
        certify_one(produced/j,('projected current saturation',label,js,rs))
        E=(2*pi/15)*norm2(Sh)*r**5
        floor=(pi/15)*(j*j/((1-g*g)**2))*r**5
        certify_one(E/floor,('projected harmonic energy floor saturation',label,js,rs))
        rows.append({'coordinate':label,'g':str(g),'M_Frobenius_squared':str(M2),'target_harmonic_shape_current':js,'r':rs,'minimum_harmonic_energy':str(E),'closed_projected_floor':str(floor)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
  'For every physical Gram coordinate g, the shape-generated STF tensor obeys |M_g|^2=2(1-g^2)^2.  Combining this identity with the exact harmonic Hodge energy floor gives E(B_r)>=(pi/15) j_h^2 r^5/(1-g^2)^2 for a harmonic shape current j_h=M_g:S_h.  Equality is attained by the pure linear harmonic strain S_h parallel to M_g, so this is the exact minimum harmonic occupancy for that coordinate current.'
 ),'rows':rows
},indent=2,allow_nan=False))
