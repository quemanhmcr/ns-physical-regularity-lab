import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS
pi=arb.pi(); sqrt2=arb(2).sqrt(); sqrt5=arb(5).sqrt()

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x > 1-t and x < 1+t):
        raise AssertionError((label,'exact ratio not tightly certified around one',x))

def certify_near_one(x,label,tol='1e-20'):
    # Asymptotic ratios are not identities: require the whole Arb interval to lie
    # inside the requested neighborhood, but do not require it to contain exactly 1.
    t=arb(tol)
    if not (x > 1-t and x < 1+t):
        raise AssertionError((label,'asymptotic ratio outside tolerance',x,t))

etas=['1e-6','0.1','1','10','1e6']
epss=['1e-6','1','1e6']
Ls=['0.5','2','10']
ks=['0','1','10','1e6','1e30']
rows=[]; asym=[]
for etas_ in etas:
  eta=arb(etas_); C=1+eta+3*eta*eta/8
  for epss_ in epss:
    eps=arb(epss_)
    for Ls_ in Ls:
      L=arb(Ls_)
      angular_reserve=((L*L+5)/(L*L)).log()/2
      if not (angular_reserve>0): raise AssertionError(('finite angular reserve not positive',etas_,epss_,Ls_))
      for ks_ in ks:
        k=arb(ks_); Q=L+eps*k*pi*C; R2=Q*Q+5
        T=-Q/(sqrt2*R2.sqrt())
        alpha=-2/R2.sqrt(); beta=1/(sqrt2*R2.sqrt())
        G=-alpha*beta
        Kba=T*alpha; Kab=-T*beta; P=Kba*Kab
        if not (Kba>0 and Kab>0 and G>0 and P>0): raise AssertionError(('strobe positive cycle lost',etas_,epss_,Ls_,ks_))
        Tlog=5*eps/(Q*R2)
        Glog=-2*eps*Q/R2
        Plog=2*eps*(5-Q*Q)/(Q*R2)
        if not (Tlog>0): raise AssertionError(('angular surplus not positive at strobe',etas_,epss_,Ls_,ks_,Tlog))
        if not (Plog-(2*Tlog+Glog)).contains(0): raise AssertionError(('strobe P rate split',etas_,epss_,Ls_,ks_))
        memory=Q-L
        expected_memory=eps*k*pi*C
        if ks_!='0': certify_one(memory/expected_memory,('renewal memory ledger',etas_,epss_,Ls_,ks_))
        long_deficit=((Q*Q+5)/(L*L+5)).log()
        rows.append({'eta':etas_,'eps':epss_,'L':Ls_,'period_index_k':ks_,'Q_memory_coordinate':str(Q),'deposited_Q_minus_L':str(memory),'T':str(T),'G_longitudinal':str(G),'P_cycle_product':str(P),'K_b_to_a':str(Kba),'K_a_to_b':str(Kab),'T_log_rate_positive_surplus':str(Tlog),'G_log_rate':str(Glog),'P_log_rate':str(Plog),'P_decreasing_after_Q_gt_sqrt5':bool(Q>sqrt5),'finite_total_angular_log_reserve':str(angular_reserve),'accumulated_longitudinal_log_deficit':str(long_deficit)})

      # At k=1e30 all tested cases are deep in the exact large-Q regime.
      k=arb('1e30'); Q=L+eps*k*pi*C; R2=Q*Q+5
      Tmag=Q/(sqrt2*R2.sqrt()); G=sqrt2/R2; P=Q*Q/(sqrt2*R2*R2)
      Kba=sqrt2*Q/R2; Kab=Q/(2*R2)
      certify_near_one(Tmag*sqrt2,('T -> 1/sqrt2',etas_,epss_,Ls_),tol='1e-20')
      certify_near_one(G*Q*Q/sqrt2,('G ~ sqrt2/Q^2',etas_,epss_,Ls_),tol='1e-20')
      certify_near_one(P*sqrt2*Q*Q,('P ~ 1/(sqrt2 Q^2)',etas_,epss_,Ls_),tol='1e-20')
      certify_near_one(Kba*Q/sqrt2,('Kba ~ sqrt2/Q',etas_,epss_,Ls_),tol='1e-20')
      certify_near_one(Kab*2*Q,('Kab ~ 1/(2Q)',etas_,epss_,Ls_),tol='1e-20')
      asym.append({'eta':etas_,'eps':epss_,'L':Ls_,'Q_at_1e30_periods':str(Q),'Tmag_times_sqrt2':str(Tmag*sqrt2),'G_Q2_over_sqrt2':str(G*Q*Q/sqrt2),'P_sqrt2_Q2':str(P*sqrt2*Q*Q),'Kba_Q_over_sqrt2':str(Kba*Q/sqrt2),'Kab_2Q':str(Kab*2*Q),'finite_total_angular_log_reserve':str(angular_reserve)})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),'asymptotic_cases':len(asym),
 'interpretation':(
   'At every periodic return of the exact quadratic NS family, d log|T|/dt is strictly positive: angular renewal surplus recurs indefinitely. '
   'Nevertheless its total logarithmic angular gain is bounded, while the exact deposited bridge-memory coordinate Q grows linearly with period count. '
   'That memory drives the longitudinal gate G like Q^-2 and the full two-edge product P like 1/(sqrt(2) Q^2), so both directed transactions die as 1/Q even though T approaches the nonzero limit 1/sqrt(2).'
 ),
 'rows':rows,'asymptotic':asym,
},indent=2,allow_nan=False))
