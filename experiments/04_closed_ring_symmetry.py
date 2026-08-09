import json, os
import mpmath as mp

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS < 160:
    raise SystemExit('ARB_PREC_BITS must be at least 160')
mp.mp.dps=max(70, int(BITS*0.30103)+30)

def vel_on_ring(R, Gamma, eps, theta):
    ct, st = mp.cos(theta), mp.sin(theta)
    xt=(R*ct,R*st,mp.mpf('0'))
    def comp(phi, idx):
        cp, sp = mp.cos(phi), mp.sin(phi)
        xp=(R*cp,R*sp,mp.mpf('0'))
        dx=(-R*sp,R*cp,mp.mpf('0'))
        r=(xt[0]-xp[0],xt[1]-xp[1],xt[2]-xp[2])
        cross=(dx[1]*r[2]-dx[2]*r[1], dx[2]*r[0]-dx[0]*r[2], dx[0]*r[1]-dx[1]*r[0])
        den=(r[0]**2+r[1]**2+r[2]**2+eps**2)**mp.mpf('1.5')
        return cross[idx]/den
    vals=[]
    twopi=2*mp.pi
    th=theta % twopi
    cuts=[mp.mpf('0')]
    for c in [th-mp.pi/2, th, th+mp.pi/2]:
        cc=c%twopi
        if cc > 0 and cc < twopi:
            cuts.append(cc)
    cuts.append(twopi)
    cuts=sorted(set(cuts))
    for idx in range(3):
        val=mp.mpf('0')
        for a,b in zip(cuts[:-1],cuts[1:]):
            if b>a:
                val += mp.quad(lambda ph: comp(ph,idx), [a,b])
        vals.append(Gamma/(4*mp.pi)*val)
    return vals

R_values=[mp.mpf('1e-3'),mp.mpf('1'),mp.mpf('1e3')]
Gamma_values=[mp.mpf('1e-6'),mp.mpf('1'),mp.mpf('1e6')]
eps_ratios=[mp.mpf('0.05'),mp.mpf('0.2'),mp.mpf('0.5')]
thetas=[mp.mpf('0'),mp.mpf('0.37'),mp.mpf('1.11'),mp.mpf('2.29'),mp.mpf('4.73')]
rel_tol=mp.mpf(10) ** (-(mp.mp.dps//3))
rows=[]
for R in R_values:
  for G in Gamma_values:
    for er in eps_ratios:
      eps=er*R
      velocities=[vel_on_ring(R,G,eps,th) for th in thetas]
      uz_ref=velocities[0][2]
      scale=max(mp.mpf('1'), abs(uz_ref))
      max_trans=max(max(abs(v[0]),abs(v[1])) for v in velocities)
      max_uz_dev=max(abs(v[2]-uz_ref) for v in velocities)
      if max_trans > rel_tol*scale:
        raise AssertionError(f'non-axial self velocity exceeds tolerance R={R},G={G},eps/R={er}: {max_trans}')
      if max_uz_dev > rel_tol*scale:
        raise AssertionError(f'axial speed not uniform around ring R={R},G={G},eps/R={er}: {max_uz_dev}')
      rows.append({'R':mp.nstr(R,12),'Gamma':mp.nstr(G,12),'eps_over_R':mp.nstr(er,12),'u_z':mp.nstr(uz_ref,30),'max_transverse_velocity':mp.nstr(max_trans,8),'max_axial_variation':mp.nstr(max_uz_dev,8)})

print(json.dumps({
  'precision_bits_requested':BITS,'mpmath_dps':mp.mp.dps,'cases':len(rows),'status':'PASS',
  'interpretation':'For an isotropically regularized circular closed vortex filament, self-induced velocity on the loop is axial and uniform, so the loop translates with zero tangential self-stretching. Strong stretching therefore requires symmetry breaking and/or nonlocal interacting geometry.',
  'rows':rows
}, indent=2))
