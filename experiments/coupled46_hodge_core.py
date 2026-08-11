import degree6_hodge_servo_core as C
from flint import arb

z=C.z; o=C.o

def degree4_basis(X,r2):
    Vb=[]; U3=[]; U5=[]; labels=[]; sectors=[]
    # Complete degree-four null poloidal sectors, carried by tangent degree-five toroidal velocities.
    for l in (1,3,5):
        fac=C.r2pow((5-l)//2,r2)
        for j,Hq in enumerate(C.harmonic_basis(l)):
            H=C.toarb(Hq); g=tuple(C.pder(H,i) for i in range(3)); T=C.cross(X,g)
            Uh=tuple(C.pmul(fac,q) for q in T); V=C.curl(Uh)
            Vb.append(V); U3.append(({}, {}, {})); U5.append(Uh); labels.append(f'P{l}_{j}'); sectors.append(('poloidal',l))
    # Complete toroidal l=4 sector with its exact tangent Hodge lift U3+U5.
    l=4; qpow=0
    A=arb(7)/22; B=-arb(4)/11; D=-arb(5)/22
    for j,Hq in enumerate(C.harmonic_basis(4)):
        H=C.toarb(Hq); g=tuple(C.pder(H,i) for i in range(3)); V=C.cross(X,g)
        low=C.vscale(D,g)
        high=tuple(C.padd(C.pscale(A,C.pmul(r2,g[i])),C.pscale(B,C.pmul(H,X[i]))) for i in range(3))
        Vb.append(V); U3.append(low); U5.append(high); labels.append(f'T4_{j}'); sectors.append(('toroidal',4))
    return Vb,U3,U5,labels,sectors

def prepare():
    X,r2,u1,u3,omega=C.setup_base()
    # Base degree-four null source.
    N4=C.sharp_split(C.bracket(omega,u3),4,X,r2)[1]
    V4b,U34,U54,L4,S4=degree4_basis(X,r2)
    K44=[]
    for V,U3 in zip(V4b,U34):
        K44.append(C.sharp_split(C.vadd(C.bracket(V,u1),C.bracket(omega,U3)),4,X,r2)[1])
    cols44=[C.flatten(v,4) for v in K44]; sel44,piv44=C.independent(cols44)
    if len(sel44)!=30: raise AssertionError(('K44 rank',len(sel44)))
    A44=[[cols44[j][piv44[i]] for j in range(30)] for i in range(30)]

    V6b,Ulow6,Uhigh6,L6,S6=C.degree6_basis(X,r2)
    K66=[C.degree6_diagonal_operator(V,low,sec,u1,omega,X,r2) for V,low,sec in zip(V6b,Ulow6,S6)]
    cols66=[C.flatten(v,6) for v in K66]; sel66,piv66=C.independent(cols66)
    if len(sel66)!=58: raise AssertionError(('K66 rank',len(sel66)))
    A66=[[cols66[j][piv66[i]] for j in range(58)] for i in range(58)]
    t4idx=[i for i,s in enumerate(S6) if s==('toroidal',4)]
    silent=[i for i,s in enumerate(S6) if s!=('toroidal',4)]
    if len(t4idx)!=9 or len(silent)!=49: raise AssertionError(('split',len(t4idx),len(silent)))
    K46=[C.degree6_lower_backreaction(low,sec,omega,X,r2) for low,sec in zip(Ulow6,S6)]
    return dict(X=X,r2=r2,u1=u1,u3=u3,omega=omega,N4=N4,
                V4b=V4b,U34=U34,U54=U54,L4=L4,S4=S4,K44=K44,cols44=cols44,piv44=piv44,A44=A44,
                V6b=V6b,Ulow6=Ulow6,Uhigh6=Uhigh6,L6=L6,S6=S6,K66=K66,cols66=cols66,piv66=piv66,A66=A66,
                t4idx=t4idx,silent=silent,K46=K46)

def solve44_field(st,target):
    b=C.flatten(target,4); x=C.solve(st['A44'],[b[p] for p in st['piv44']])
    V=C.combine(x,st['V4b']); U3=C.combine(x,st['U34']); U5=C.combine(x,st['U54'])
    res=C.vadd(C.combine(x,st['K44']),C.vscale(-1,target))
    return x,V,U3,U5,res

def solve66_field(st,target):
    b=C.flatten(target,6); x=C.solve(st['A66'],[b[p] for p in st['piv66']])
    V=C.combine(x,st['V6b']); res=C.vadd(C.combine(x,st['K66']),C.vscale(-1,target))
    return x,V,res

def v4_for_y(st,y):
    # y lives in the physical degree-six T4 feedback sector.
    V6y=C.combine(y,[st['V6b'][i] for i in st['t4idx']])
    U3y=C.combine(y,[st['Ulow6'][i] for i in st['t4idx']])
    B4=C.combine(y,[st['K46'][i] for i in st['t4idx']])
    rhs=C.vscale(-1,C.vadd(st['N4'],B4))
    c4,V4,U34,U54,res=solve44_field(st,rhs)
    return dict(V6y=V6y,U3y=U3y,B4=B4,c4=c4,V4=V4,U34=U34,U54=U54,res4=res)

def nonlinear_degree6_source(st,y):
    q=v4_for_y(st,y); V4=q['V4']; U34=q['U34']; U54=q['U54']; U3y=q['U3y']
    R=C.vadd(C.bracket(V4,st['u3']),C.bracket(st['omega'],U54))
    R=C.vadd(R,C.bracket(V4,U34))
    R=C.vadd(R,C.bracket(V4,U3y))
    return C.sharp_split(R,6,st['X'],st['r2'])[1],q

def feedback_map(st,y):
    R,q=nonlinear_degree6_source(st,y)
    coeff,V6,res=solve66_field(st,C.vscale(-1,R))
    phi=[coeff[i] for i in st['t4idx']]
    F=[phi[j]-y[j] for j in range(9)]
    W=C.combine([coeff[i] for i in st['silent']],[st['V6b'][i] for i in st['silent']])
    return dict(phi=phi,F=F,coeff6=coeff,V6=V6,W=W,R6=R,res6=res,**q)

def feedback_jacobian(st,y):
    """Exact directional Jacobian of F(y)=Pi_T4 K66^{-1}[-R6(V4(y),y)]-y."""
    q=v4_for_y(st,y); V4=q['V4']; U34=q['U34']; U3y=q['U3y']
    cols=[]
    for a,idx in enumerate(st['t4idx']):
        # dV4/dy_a is fixed by exact lower-level compensation K44 dV4 = -K46 e_a.
        _,dV4,dU3,dU5,_=solve44_field(st,C.vscale(-1,st['K46'][idx]))
        dU3y=st['Ulow6'][idx]
        dR=C.vadd(C.bracket(dV4,st['u3']),C.bracket(st['omega'],dU5))
        dR=C.vadd(dR,C.bracket(dV4,q['U34']))
        dR=C.vadd(dR,C.bracket(V4,dU3))
        dR=C.vadd(dR,C.bracket(dV4,U3y))
        dR=C.vadd(dR,C.bracket(V4,dU3y))
        dR=C.sharp_split(dR,6,st['X'],st['r2'])[1]
        dc,_,_=solve66_field(st,C.vscale(-1,dR))
        cols.append([dc[j]-(o if i==a else z) for i,j in enumerate(st['t4idx'])])
    # Return row-major 9x9 Jacobian.
    return [[cols[j][i] for j in range(9)] for i in range(9)]

def pswap_xy(P):
    return {(e[1],e[0],e[2]):v for e,v in P.items()}

def vpolar_swap_xy(V):
    # Polar-vector action of the improper orthogonal reflection x<->y.
    return (pswap_xy(V[1]),pswap_xy(V[0]),pswap_xy(V[2]))

def vaxial_swap_xy(V):
    # Vorticity is axial: det(R) R with det(R)=-1 for x<->y.
    return C.vscale(-1,vpolar_swap_xy(V))

def feedback_symmetry_basis(st):
    """Return the physical x<->y fixed/anti-fixed decomposition of the 9D degree-six T4 feedback sector."""
    Y=[st['V6b'][i] for i in st['t4idx']]
    cols=[C.flatten(v,6) for v in Y]; sel,piv=C.independent(cols)
    if len(sel)!=9: raise AssertionError(('T4 basis rank',len(sel)))
    A=[[cols[j][piv[i]] for j in range(9)] for i in range(9)]
    swapcols=[]
    for V in Y:
        b=C.flatten(vaxial_swap_xy(V),6); c=C.solve(A,[b[p] for p in piv]); swapcols.append(c)
    sym=[]; anti=[]
    for j,c in enumerate(swapcols):
        ej=[z]*9; ej[j]=o
        sym.append([ej[i]+c[i] for i in range(9)])
        anti.append([ej[i]-c[i] for i in range(9)])
    ss,_=C.independent(sym); aa,_=C.independent(anti)
    S=[sym[i] for i in ss]; Aanti=[anti[i] for i in aa]
    if len(S)!=5 or len(Aanti)!=4: raise AssertionError(('swap fixed dimensions',len(S),len(Aanti)))
    _,spiv=C.independent(S)
    Smat=[[S[j][spiv[i]] for j in range(5)] for i in range(5)]
    return dict(Y=Y,swapcols=swapcols,S=S,Aanti=Aanti,spiv=spiv,Smat=Smat)

def y_from_sym(a,sym):
    return [sum((a[k]*sym['S'][k][i] for k in range(5)),z) for i in range(9)]

def sym_coords(y,sym):
    a=C.solve(sym['Smat'],[y[p] for p in sym['spiv']])
    recon=y_from_sym(a,sym)
    return a,recon

def reduced_feedback(st,sym,a):
    y=y_from_sym(a,sym); fb=feedback_map(st,y); g,recon=sym_coords(fb['F'],sym)
    return g,fb,recon

def reduced_jacobian(st,sym,a):
    y=y_from_sym(a,sym); J9=feedback_jacobian(st,y)
    cols=[]
    for k in range(5):
        dy=sym['S'][k]
        df=[sum((J9[i][j]*dy[j] for j in range(9)),z) for i in range(9)]
        c,_=sym_coords(df,sym); cols.append(c)
    return [[cols[j][i] for j in range(5)] for i in range(5)]

def matvec(A,x):
    return [sum((A[i][j]*x[j] for j in range(len(x))),z) for i in range(len(A))]

def matmul(A,B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))),z) for j in range(len(B[0]))] for i in range(len(A))]

def eye(n):
    return [[o if i==j else z for j in range(n)] for i in range(n)]
