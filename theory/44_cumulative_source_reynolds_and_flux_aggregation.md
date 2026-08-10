# Cumulative source-Reynolds dichotomy: renewal exposure or circulation-flux aggregation

## THINK — one high-Re event is not enough

The localized spectral gap gives a per-production-time persistence parameter

`Re_source=s L^2/nu`.

But a finite-time vorticity blow-up requires infinite cumulative positive production, so the relevant clock is the accumulated number of strain times

`dN=s dt`, `N->infinity`.

The localized viscous exposure along that history is

`dXi=(nu/L^2)dt=dN/Re_source`.

Thus

`Xi=integral dN/Re_source`.

A source can have `Re_source>>1` at every instant and still accumulate infinite viscous exposure if it executes infinitely many production times.

## Exact cumulative dichotomy

For any fixed threshold `R_*`,

`integral_{Re_source<=R_*} dN`

`<= R_* Xi`.

Therefore if the total strain count diverges while `Xi<infinity`, only finite strain count can occur in every bounded-Re sector.  The source must visit arbitrarily large Reynolds number:

`limsup Re_source=infinity`.

Since

`Gamma_Q=sL^2=nu Re_source`,

finite cumulative exposure forces unbounded circulation-dimensional transaction along a subsequence.

This is stronger than saying the branch is “high Re”: it says a genuinely frozen escape through an infinite cascade must become **more and more high Re**.

## Physical flux of a localized extremal source

For the sharp extremal ray `E=diag(2,-1,-1)` and localized radial profile

`q=q0 x^2(1-x^p)`, `x=r/L`,

the Hodge source strain at `R=L` is

`s=C_p q0`,

`C_p=5p(p+9)/[14(p+2)(p+7)]`.

The exact Clebsch flux cell has

`Delta A=5p q0/[6(p+2)]`,

`Delta B=(3/2)L^2`,

so

`Phi_ext=Delta A Delta B`

`=[7(p+7)/(2(p+9))] sL^2`

`=[7(p+7)/(2(p+9))] nu Re_source`.

Hence the unbounded-Re escape is exactly an **unbounded physical vorticity-flux aggregation** in the minimum productive carrier, not merely a large dimensionless number.

## ATTACK — scalar budgets do not kill this branch

A scale ledger demonstrates why another scalar estimate is unlikely to finish the proof.  Parameterize an inward cascade by strain count `N` and choose

`Re_source=e^(aN)`,

`L=e^(-bN)`,

with `b>2a>0` and `nu=1` for units.

Then

`s=Re/L^2=e^[(a+2b)N]`,

so infinite `N` fits into finite physical time.  Moreover:

- remaining localized viscous exposure is `~e^(-aN)` and finite;
- the natural strain occupancy scale `s^2L^5=e^[(2a-b)N]` decays;
- the modeled sharp-enstrophy dissipation tail is `~e^[(a-b)N]` and finite;
- but `Phi_ext~Re=e^(aN)` diverges.

This is **not an exact Navier-Stokes solution**.  It is an escape-route calibration showing that all presently identified scalar energy/dissipation scalings can coexist with a frozen high-Re inward cascade if ancestry flux is allowed to aggregate without geometric restriction.

Therefore the missing obstruction is not “a stronger scalar toll.”  It is a theorem about how material ancestry can be recruited, packed, folded, or merged into a shrinking productive source.

## PROMOTE / KILL

PROMOTE:
1. The correct cumulative persistence clock is `Xi=int dN/Re_source`.
2. Infinite strain count with finite `Xi` forces `Re_source` and `Gamma_Q` unbounded along a subsequence.
3. On a localized sharp carrier this is exactly unbounded physical circulation-flux aggregation.
4. The moving-source ancestry ledger identifies relative material recruitment as the high-Re supply channel complementary to viscous renewal.

KILL:
1. “Re_source is large at each event” as sufficient evidence that ancestry stays frozen through an infinite cascade.
2. A proof strategy relying only on the currently known scalar energy/enstrophy/dissipation costs in the high-Re branch.
3. A fixed circulation toll per scale; the actual issue is aggregation of flux ancestry into a shrinking source.

## Next frontier — material aggregation geometry

The remaining high-Re question is now concrete:

**Can a shrinking productive source aggregate unbounded approximately frozen circulation flux through its moving boundary while staying near the sharp productive carrier and avoiding unbounded transaction-null/folding geometry?**

That is the next closed-network/material-packing microscope.
