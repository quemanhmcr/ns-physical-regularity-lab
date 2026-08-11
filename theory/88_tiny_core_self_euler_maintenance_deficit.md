# Tiny screened core cannot maintain itself

The compact physical-space collar of module 159 is cheap in enstrophy and palinstrophy.  Could it use its own induced velocity to replenish itself against viscosity?

Let

`omega_epsilon=a_epsilon(r) x cross grad H_l`,

with `C_l[a_epsilon]=1` and core radius `epsilon`.

Inside the core, the exact tangent Hodge velocity is obtained from the radial div-curl system.  In scaled coordinates `x=epsilon y`, its structure is

`omega_epsilon = epsilon^(l-2) Omega_epsilon(y)`,

`u_epsilon = epsilon^(l-1) U_epsilon(y)`.

The scaled shapes remain finite because the normalization amplitude `A_epsilon` tends to a finite limit.

Therefore the actual self-Euler vorticity source scales as

`(omega.grad)u-(u.grad)omega = epsilon^(2l-4) G_epsilon(y)`,

whereas

`Delta omega = epsilon^(l-4) D_epsilon(y)`.

The full source-norm ratio is consequently

`||Euler_self|| / ||nu Delta omega|| = O(epsilon^l/nu)`.

More importantly, the same scaling survives the actual screened source ledger.  Project the self-Euler source back into the original angular `T_l` channel.  Its radial coefficient carries the prefactor `epsilon^(l-4)`, and the screened radial integral contributes `epsilon^2`, so

`J_0^self = O(epsilon^(l-2))`.

But the viscous drift is exact:

`M_1=C_l[D_l a_epsilon]=(l+1)A_epsilon epsilon^-2`.

Thus

`J_0^self/M_1 = O(epsilon^l)`.

For every `l>=2`, this tends to zero.

So the tiny core escape is not a closed microscopic machine.  The same structure that makes its static enstrophy and palinstrophy cheap also makes its self-induced Euler turnover too weak to oppose viscosity at sufficiently small scale.

Continuous maintenance must therefore be supplied by another actor:

- a strain field whose rate grows like `nu epsilon^-2` (modules 158-160), or
- genuinely non-affine interactions with other same-scale vorticity/velocity structures.

The second option is now the important frontier because it necessarily introduces additional material ancestry.  A single isolated screened core cannot recycle itself indefinitely.
