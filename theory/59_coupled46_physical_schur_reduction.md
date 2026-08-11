# Coupled degree-four/degree-six Hodge feedback reduction

## THINK — do not force a sequential hierarchy on a triangular physical system

The degree-six diagonal servo is unique, but its toroidal `l=4` component carries a harmonic degree-three velocity companion and reopens degree four.  Therefore the sequential picture

`solve degree 4 -> solve degree 6`

is not the native organization of the Hodge/Euler response.

The correct question is the simultaneous degree-four/degree-six compatibility problem.

## Physical split of the degree-six control space

The complete degree-six null Hodge space is

`N_6 = P1 + P3 + P5 + P7 + T4 + T6`,

of dimension `58`.

Only the `T4` sector carries a harmonic velocity companion of degree three.  Therefore only `T4` can act on the base degree-two vorticity and produce a degree-four response.

Thus nature supplies the intrinsic decomposition

`N_6 = Y + W`,

where

`Y=T4`, `dim Y=9`,

and

`W=P1+P3+P5+P7+T6`, `dim W=49`.

The split is not chosen by an algebraic convenience: it is determined by which Hodge velocity companions can physically reach the lower response level.

## Exact elimination of degree four

Let `y in Y` be the degree-six feedback component.  Its lower Hodge companion creates a degree-four response `K46 y`.

Because the complete degree-four Hodge operator `K44` is invertible, degree-four cancellation uniquely determines

`V4(y) = -K44^{-1}(N4 + K46 y)`.

So degree four is not an independent unknown once the physical feedback `y` is specified.

## Exact elimination of the lower-silent degree-six sector

For fixed `y`, the degree-six null source produced by `V4(y)` and the lower companion of `y` is

`R6(V4(y),y)`.

The complete degree-six diagonal operator `K66` is invertible, so

`V6 = K66^{-1}[-R6(V4(y),y)]`

is uniquely determined.

Write its physical split as

`V6 = Pi_Y V6 + Pi_W V6`.

The 49-dimensional `W` part is automatically the unique lower-silent response required by the degree-six equation.  Coupled consistency requires only

`y = Pi_Y V6`.

Therefore the full simultaneous problem reduces exactly to

`y = Pi_T4 K66^{-1}[-R6(V4(y),y)]`.

This is a nine-dimensional physical feedback law.

## Why this reduction matters

The apparent `30+58` response system was partly observer complexity.  Once Hodge companion degrees are respected, only nine feedback directions communicate from degree six back to degree four.  All other amplitudes are slaved linearly.

The next question is not whether an 88-dimensional nonlinear solver finds a root.  It is whether this intrinsic nine-channel feedback map has a fixed point.

## ATTACK criteria

Module 131 must remotely certify:

1. `rank K46|_{T4}=9`;
2. every vector in the 49-dimensional complement is structurally lower-silent;
3. `rank K66|_W=49`;
4. exact degree-four elimination works for the sequential-feedback calibration;
5. the old sequential solution is generally not already a fixed point of the coupled feedback map.

Only after these facts survive 160/256/512 precision is it legitimate to solve the nonlinear nine-dimensional law.

## Do not promote prematurely

This reduction by itself does not prove existence or nonexistence of a coupled servo.  It only identifies the correct physical coordinates in which that question should be asked.
