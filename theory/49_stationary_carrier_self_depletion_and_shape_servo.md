# Stationary productive amplification needs a shape servo and a sacrificial carrier collar

## THINK — can the critical winding donor generate the tensor that keeps itself productive?

The winding escape is scalar unless the full pair Gram geometry can remain productive.  At the capacity pair geometry, prescribe endpoint vorticity gains `s_a,s_b>0` while demanding

`alphadot=betadot=gammadot=0`.

The pair-generated STF frame determines the common strain uniquely.

## Exact asymmetric stationary-lock tensor

The tensor is

`S_lock = [[s_a,0,sqrt(2)(2s_a+s_b)/4],`

`          [0,s_b,sqrt(2)(s_a+2s_b)/4],`

`          [...,...,-s_a-s_b]]`.

It satisfies the requested endpoint gains and zero Gram-shape rates.  Its invariants are

`|S_lock|^2=[13s_a^2+16s_as_b+13s_b^2]/4`,

and

`det S_lock`

`=-(s_a+s_b)(s_a^2+15s_as_b+s_b^2)/8 <0`.

So every positive two-endpoint stationary shape lock lies in the negative-determinant strain sector.
For extreme asymmetric gains, the small endpoint response must be certified directly in pair-frame coordinates.  Reconstructing `s_a<<s_b` by contracting the assembled Cartesian tensor subtracts order-`s_b` parent terms and creates a false precision problem even though the response-frame identity is exact.
  The bridge direction contracts at rate

`n.S_lock.n=-(s_a+s_b)`,

which is exactly the pair-cell determinant balance.

## The minimum tangent carrier does not globally self-amplify

For any STF strain `S`, use its exact self-contained tangent Hodge carrier

`u_c=(1-5r^2/(3L^2))Sx+[2/(3L^2)](x.S.x)x`,

`omega_c=-(14/(3L^2))x cross Sx`.

Direct invariant sphere moments give the shell-averaged stretching coefficient `P_shell/Z_shell`

`c_shell(r)`

`=(3/7)[det S/|S|^2][13(r/L)^2-7]`.

After volume integration,

`[integral omega_c.S_{u_c}omega_c dx]/[integral |omega_c|^2 dx]`

`=(4/3) det S/|S|^2`.

Therefore the minimum sharp carrier of every positive stationary-lock tensor has **negative total self-stretch**.

In the symmetric case `s_a=s_b=s`,

`det S=-17s^3/4`, `|S|^2=21s^2/2`,

so

`P_self/Z=-34s/63`.

In the enstrophy equation `(1/2) Zdot=P_self-...`, the self-stretch contribution to `d log Z/dt` is twice this coefficient, `-68s/63`.

The shell rate changes sign at

`r/L=sqrt(7/13)`.

For negative determinant, the inner carrier core is self-amplified while the outer collar is self-deamplified more strongly.  This is a core-gain/collar-loss stretching ledger, not a conserved transfer statement.

Thus the exact minimum carrier is not a globally self-sustaining winding engine.  Sustained critical operation must replenish/rearrange the outer productive ancestry or use additional vorticity outside the minimum carrier.

## Two production channels do not close the shape dynamics

A second tensorial attack is even sharper.  At the symmetric capacity geometry, take only the two endpoint-axis STF production rays

`A_a=a tensor a-I/3`,

`A_b=b tensor b-I/3`.

Choosing

`Q_2=(9s/4)(A_a+A_b)`

gives exactly the desired endpoint gains `s,s`.  But its Gram-shape rates are

`alphadot=2sqrt(3)s/3`,

`betadot=-2sqrt(3)s/3`,

`gammadot=-4s/3`.

So the scalar positive two-cycle immediately drifts out of its locked geometry.

The exact correction is

`Q_servo=-3s A_n`,

where

`A_n=n tensor n-I/3`.

At the capacity geometry,

`A_a:Q_servo=A_b:Q_servo=0`,

while its three shape rates are exactly the negatives of the two-ray drift.  Hence

`S_lock=Q_2+Q_servo`

in this STF response representation.

This does **not** mean nature contains three independent abstract graph edges.  It means the complete physical transaction tensor has an indispensable shape-servo component that is invisible to the two magnitude-production channels.

## PROMOTE / KILL

PROMOTE:
1. Positive stationary pair amplification universally requires a negative-determinant lock tensor.
2. The sharp self-contained carrier of that tensor has negative global self-stretch, with inner core gain and stronger outer collar depletion.
3. Matching the two endpoint production rates alone necessarily drives Gram-shape drift.
4. Shape-lock requires an additional pure shape transaction component that does no endpoint magnitude work at the capacity geometry.

KILL / DEMOTE:
1. A scalar two-node positive-cycle ODE as a closed physical model of the critical winding branch.
2. The minimum tangent carrier as a globally self-sustaining source of its own stationary amplifier.

## Next frontier — replenishment of the sacrificial collar

The critical winding escape can survive only if through-going ancestry continually replenishes the outer productive carrier/shape-servo sector while the inner core is amplified.  The remaining question is therefore not simply whether two lineages mutually stretch, but whether a closed material network can maintain the required **five-channel tensor response** and replenish its negative-self-stretch collar indefinitely.

That is a network current problem with a genuine source/sink structure in physical transaction channels.
