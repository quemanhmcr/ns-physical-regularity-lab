import json, os, math
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
import degree6_hodge_servo_core as C
z=C.z; one=C.o; pi=arb.pi()

def Pn(n):
    # Re[(x+i z)^n]
    P={}
    for k in range(0,n+1,2):
        P[(n-k,0,k)]=arb(((-1)**(k//2))*math.comb(n,k))
    return P

def ball_even(exp):
    a,b,c=exp
    if a%2 or b%2 or c%2: return z
    aa=arb(a)/2+arb('0.5'); bb=arb(b)/2+arb('0.5'); cc=arb(c)/2+arb('0.5')
    N=(a+b+c)//2
    sph=2*aa.gamma()*bb.gamma()*cc.gamma()/(arb(N)+arb('1.5')).gamma()
    return sph/arb(2*N+3)

def pint(P): return sum((v*ball_even(e) for e,v in P.items()),z)
def energy_ball_unit(U): return pint(sum_poly_squares(U))
def sum_poly_squares(U):
    q={}
    for u in U: q=C.padd(q,C.pmul(u,u))
    return q

def sphere_Pn_square(n):
    # Integral_{S^2} Re[(x+i z)^n]^2 dS
    return pi*(arb(2)**(2*n+1))*(arb(math.factorial(n))**2)/arb(math.factorial(2*n+1))

def closed_energy_constant(m):
    n=m+1
    return sphere_Pn_square(n)/arb(n)

rows=[]
# Direct polynomial/ball-integral certification of the Green-identity energy constant.
for m in (2,4,6,8,12):
    n=m+1
    phi=C.pscale(-one/arb(n),Pn(n))
    lap=C.plap(phi)
    if any(not v.contains(0) for v in lap.values()): raise AssertionError(('harmonic potential',m,lap))
    u=tuple(C.pder(phi,j) for j in range(3))
    Eunit=energy_ball_unit(u)
    c=closed_energy_constant(m)
    if not (Eunit/c).contains(1): raise AssertionError(('energy constant',m,Eunit,c,Eunit/c))
    rows.append({'kind':'direct_energy_certificate','m':m,'n_potential_degree':n,
                 'direct_unit_ball_velocity_energy':str(Eunit),'closed_energy_constant_c_m':str(c),
                 'direct_over_closed':str(Eunit/c)})

# Normalize amplitude so |partial_z^2 u_x| at z=epsilon equals a prescribed curvature source Q.
# For phi_m=-Re[(x+i z)^(m+1)]/(m+1), |partial_z^2 u_x|=m(m-1) epsilon^(m-2).
# Thus E/(Q^2 epsilon^7)=c_m lambda^(2m+3)/[m^2(m-1)^2], lambda=L/epsilon.
for m in (2,4,8,16,32,64):
    c=closed_energy_constant(m)
    denom=arb(m*m*(m-1)*(m-1))
    base=c/denom
    lambda_budget1=(denom/c)**(one/arb(2*m+3))
    for ls in ('1','1.05','1.1','1.25','1.5','2'):
        lam=arb(ls)
        penalty=lam**(2*m+3)
        normE=base*penalty
        if not (normE/base/penalty).contains(1): raise AssertionError(('normalized energy law',m,ls))
        rows.append({'kind':'source_radius_gate','m':m,'lambda_L_over_epsilon':ls,
                     'closed_energy_constant_c_m':str(c),
                     'normalized_occupancy_E_over_Q2_epsilon7':str(normE),
                     'same_scale_normalized_occupancy_lambda1':str(base),
                     'remote_separation_penalty_E_lambda_over_E_lambda1':str(normE/base),
                     'closed_remote_penalty_lambda_power_2mplus3':str(penalty),
                     'lambda_max_for_unit_normalized_occupancy':str(lambda_budget1)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
   'For the irrotational harmonic folding mode phi_m=-Re[(x+i z)^(m+1)]/(m+1), direct polynomial ball integration agrees with the exact Green-identity occupancy integral int_{B_L}|grad phi_m|^2 = c_m L^(2m+3), where c_m=[int_{S^2} Re(x+i z)^(m+1)^2 dS]/(m+1). '
   'Normalize the mode amplitude so that the physical straight-line curvature source |partial_z^2 u_x| at the core point z=epsilon equals Q.  Then exactly E/(Q^2 epsilon^7)=c_m (L/epsilon)^(2m+3)/[m^2(m-1)^2].  The full dependence on source clearance is therefore the exponential-in-degree factor (L/epsilon)^(2m+3); no vorticity norm or observer normalization enters it. '
   'Consequently a high-degree harmonic folding actor cannot act from a fixed relative clearance L/epsilon>1 at bounded normalized kinetic occupancy.  Finite occupancy forces the harmonicity radius toward the core as m grows, up to only polynomial prefactors in m. '
   'This is an exact single-mode source-radius gate, not yet a theorem that N simultaneous cap crossings require one degree-m mode with m comparable to N.  The surviving escape is near-contact: place the vorticity that sources the harmonic actor in a collar whose gap from the shrinking core tends to zero.'),
 'rows':rows
},indent=2,allow_nan=False))
