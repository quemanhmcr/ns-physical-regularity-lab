# Eigenframe symmetry explains the missing P1 channel

The generic nonsymmetric feedback attack kills the hypothesis that `P1` is universally protected.  The certified coupled branch nevertheless has zero `P1`.  The correct explanation is the stabilizer of the stationary strain.

A real symmetric strain tensor with distinct eigenvalues is invariant under the three proper pi rotations about its orthogonal eigenaxes.  For the capacity strain, one eigenaxis is `(1,-1,0)`, giving

`R1(x,y,z)=(-y,-x,-z)`.

A second eigenaxis lies in the span of `(1,1,0)` and `ez`, proportional to

`e+ +(sqrt(2)-1) ez`.

Let `R2` be the pi rotation about this axis.  The common fixed-vector space of `R1` and `R2` is `{0}`.

All constructions used to obtain the local Hodge servo—curl, divergence, the ball Hodge split, the transaction projector, Euler brackets and the fixed linear response inverses—are equivariant under proper rotations.  The branch generated from the stationary carrier should therefore inherit its eigenframe stabilizer.

Module 141 checks this directly on the physical polynomial fields `u1,u3,omega2,V4,V6,N8,N10,N12`.  If they are invariant under both rotations, their sphere-mean vorticity vectors must lie in the common fixed-vector space and hence vanish.  Since the homogeneous `P1` sector is exactly the sphere-mean-bearing poloidal sector, this gives the physical selection rule

`D2 eigenframe symmetry -> P1 = 0`.

This is a symmetry theorem for the stationary-capacity branch, not a universal conservation law for Navier-Stokes.
