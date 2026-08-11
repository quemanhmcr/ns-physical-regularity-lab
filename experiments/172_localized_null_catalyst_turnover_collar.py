import json, os
from fractions import Fraction as F
from math import comb
from flint import arb, ctx
BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi()

def A(q): return arb(q.numerator)/q.denominator

def poly_chi(p):
    # chi(s)=(1-s^2)^p, 0<=s<=1.
    return {2*k:F(((-1)**k)*comb(p,k)) for k in range(p+1)}
def pint_product(c,d,power):
    out=F(0)
    for i,a in c.items():
        for j,b in d.items(): out += a*b/F(power+i+j+1)
    return out
def deriv(c): return {k-1:F(k)*v for k,v in c.items() if k>0}
def chi_value(c,s):
    out=arb(0)
    for k,v in c.items(): out += A(v)*(s**k)
    return out

rows=[]
for p in (2,3,4,8,16):
    c=poly_chi(p); dc=deriv(c)
    IE=pint_product(c,c,6)
    ID=pint_product(dc,dc,6)
    # Raw enstrophy angular formula divided by 8 pi A^2 L^5:
    # int s^4 [chi^2 +(2/5)s chi chi' +(1/15)s^2 chi'^2] ds.
    raw1=pint_product(c,c,4)
    cross=pint_product(c,dc,5)*F(2,5)
    grad=pint_product(dc,dc,6)*F(1,15)
    raw=raw1+cross+grad
    reduced=ID*F(1,15)
    if raw!=reduced: raise AssertionError(('turnover identity',p,raw,reduced))
    # E=(4pi/15)A^2 L^7 IE, Z=(8pi/15)A^2 L^5 ID.
    clock_coeff=F(1,2)*IE/ID # E/(nu Z) = clock_coeff L^2/nu.
    # radial cap flux is proportional to s^3 chi(s), maximum at s^2=3/(2p+3).
    speak=(arb(3)/arb(2*p+3)).sqrt(); fpeak=(speak**3)*chi_value(c,speak)
    rows.append({
        'taper_power_p':p,'chi_form':'(1-s^2)^p on [0,1], zero outside',
        'dimensionless_energy_integral_int_s6_chi2':str(A(IE)),
        'dimensionless_turnover_integral_int_s6_chiprime2':str(A(ID)),
        'raw_enstrophy_integral_before_IBP':str(A(raw)),
        'reduced_enstrophy_integral_one_over15_ID':str(A(reduced)),
        'energy_over_nu_enstrophy_clock_coefficient_times_L2_over_nu':str(A(clock_coeff)),
        'peak_radial_flux_radius_over_L':str(speak),
        'peak_radial_cap_flux_shape_s3_chi':str(fpeak),
    })

# Sharp weighted Cauchy infimum for any turnover from chi(R)=1 to chi(infty)=0:
# int_R^inf r^6 chi'^2 dr >= 1/int_R^inf r^-6 dr =5 R^5.
# Equality profile in the relaxed H1 class is chi=(R/r)^5 outside R.
R=arb(1)
sharp_ID=arb(5)*R**5
# For B norm^2=6: Z >= (8pi/15) A^2 *5 R^5=(8pi/3)A^2 R^5.
sharp_Z_coeff=arb(8)*pi/3
# Relaxed equality energy: inside chi=1 plus r^-5 tail, int r^6 chi^2 dr=R^7(1/7+1/3)=10R^7/21.
sharp_IE=arb(10)/21
sharp_clock=(arb(1)/2)*sharp_IE/sharp_ID # R=1

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'sharp_turnover_weighted_derivative_infimum_coefficient_5R5':str(sharp_ID),
 'Bdiag_enstrophy_infimum_coefficient_Z_over_A2R5':str(sharp_Z_coeff),
 'relaxed_infimum_energy_integral_coefficient_10over21_R7':str(sharp_IE),
 'relaxed_infimum_energy_dissipation_clock_coefficient_R2_over_nu':str(sharp_clock),
 'interpretation':(
   'Localize the validated linear null catalyst by multiplying its exact tangent velocity u0 by a radial cutoff chi(r/L), then define omega=curl(chi u0).  This preserves divergence-free vorticity automatically.  For B=diag(2,-1,-1), exact angular integration gives E=(4 pi A^2/15) int r^6 chi^2 dr. '
   'The full enstrophy initially contains bulk chi^2, cross chi chi_prime, and derivative terms, but integration by parts cancels the bulk exactly and leaves Z=(8 pi A^2/15) int r^6 chi_prime^2 dr.  Thus the unlocalized linear germ is viscosity-null only because it has no turnover; every finite-energy realization puts its viscous exposure in the spatial closure collar. '
   'The radial vorticity flux through centered spheres is multiplied by chi, so it rises as r^3 chi(r/L) and returns to zero at the outer boundary: the same collar is also where the through-going null ancestry turns over. '
   'For fixed taper shape the kinetic-energy dissipation clock E/(nu Z) is an order-one constant times L^2/nu.  A broad halo therefore can be reused for order (L/epsilon)^2 core viscous clocks; localization restores viscosity but does not by itself kill the broad-reservoir escape. '
   'A weighted Cauchy calculation gives the sharp relaxed turnover infimum int_R^inf r^6 chi_prime^2 >=5R^5.  This is a physical closure-collar derivative floor, not yet an irreversible per-event cost.'),
 'rows':rows
},indent=2,allow_nan=False))
