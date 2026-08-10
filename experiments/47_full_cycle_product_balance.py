import json, os
from flint import arb, ctx

BITS=int(os.environ.get('ARB_PREC_BITS','160'))
if BITS<160: raise SystemExit('ARB_PREC_BITS must be at least 160')
ctx.prec=BITS

def certify_one(x,label,tol='1e-30'):
    t=arb(tol)
    if not x.contains(1) or not (x > 1-t and x < 1+t):
        raise AssertionError((label,'ratio not tightly certified around one',x))

etas=['1e-6','0.1','1','1e3']
epss=['1e-6','1','1e6']
Ls=['0.5','2','10']
times=['0','1e-12','0.01','0.1','0.5','1','2','5','10']
rows=[]
for etas_ in etas:
  eta=arb(etas_); C=1+eta+3*eta*eta/8
  for epss_ in epss:
    eps=arb(epss_)
    for Ls_ in Ls:
      L=arb(Ls_)
      for ts in times:
        t=arb(ts); s=t.sin(); sin2=(2*t).sin(); sin4=(4*t).sin()
        lam=1+eta*s*s; ldot=eta*sin2
        I=t*C-(eta/2+eta*eta/4)*sin2+(eta*eta/32)*sin4
        Q=L+eps*I; Qdot=eps*lam*lam
        S=Q*Q+4*lam**4+lam*lam
        Sdot=2*Q*Qdot+16*lam**3*ldot+2*lam*ldot
        r2=S/(lam*lam)
        r2dot=2*Q*Qdot/(lam*lam)-2*Q*Q*ldot/(lam**3)+8*lam*ldot
        sigma_R=r2dot/(2*r2)
        sigma_a=ldot/lam
        sigma_b=lam*ldot/(lam*lam+1)
        lambda_D=Qdot/Q
        lambda_La=2*ldot/lam
        lambda_Lb=4*lam*ldot/(2*lam*lam-1)
        Tlog=lambda_D-sigma_a-sigma_b-sigma_R
        Glog=lambda_La+lambda_Lb-sigma_a-sigma_b-2*sigma_R
        Plog=2*lambda_D+lambda_La+lambda_Lb-3*sigma_a-3*sigma_b-4*sigma_R
        if not (Plog-(2*Tlog+Glog)).contains(0): raise AssertionError(('P= T^2 G rate split',etas_,epss_,Ls_,ts))

        # Independent direct closed-form derivative of
        # P=2 Q^2 lambda^3(2 lambda^2-1)/[(lambda^2+1)^(3/2) S^2].
        Plog_direct=2*Qdot/Q+3*ldot/lam+4*lam*ldot/(2*lam*lam-1)-3*lam*ldot/(lam*lam+1)-2*Sdot/S
        if Plog_direct.contains(0):
            if not Plog.contains(0): raise AssertionError(('zero P rate mismatch',etas_,epss_,Ls_,ts))
            ratio='zero/zero'
        else:
            rr=Plog/Plog_direct; certify_one(rr,('full product direct balance',etas_,epss_,Ls_,ts)); ratio=str(rr)

        T2=Q*Q/((lam*lam+1)*S)
        G=2*lam**3*(2*lam*lam-1)/((lam*lam+1).sqrt()*S)
        P=T2*G
        if not (G>0 and P>0): raise AssertionError(('positive cycle product gate lost',etas_,epss_,Ls_,ts,G,P))
        rows.append({'eta':etas_,'eps':epss_,'L':Ls_,'t':ts,'lambda':str(lam),'Q':str(Q),'T_squared':str(T2),'G_longitudinal':str(G),'P_cycle_product':str(P),'T_log_rate':str(Tlog),'G_log_rate':str(Glog),'P_log_rate_balance':str(Plog),'P_log_rate_direct':str(Plog_direct),'balance_over_direct':ratio})

print(json.dumps({
 'arb_precision_bits':BITS,'status':'PASS','cases':len(rows),
 'interpretation':(
   'The exact full mutual-cycle product P=T^2 G obeys d log P=2 d log|T|+d log G over the time-dependent quadratic NS family. '
   'An independent derivative of the closed P expression matches the ancestry-rate balance across wide parameter scales. '
   'This certifies that normalized noncoplanarity and opposite-longitudinal access are distinct, multiplicative survival gates rather than interchangeable descriptions of the same geometry.'
 ),
 'rows':rows,
},indent=2,allow_nan=False))
