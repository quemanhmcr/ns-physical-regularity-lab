# NS Physical Regularity Lab

Physics-first laboratory for stress-testing a proposed Navier–Stokes regularity mechanism:

> A singularly viable transfer to smaller physical scales cannot occur for free; each essential scale transition must either incur an irreversible physical toll or lose the coherence/circulation needed to keep outrunning viscosity.

## Research discipline

- Navier–Stokes itself is the evolution law. No artificial evolution is introduced.
- Start from physical mechanisms (material deformation, vorticity stretching, incompressibility, circulation, pressure response, viscosity) before abstract norm closure.
- Numerical computation is **never run locally**. Linux is used only to edit, version, and submit code.
- All experiments execute on GitHub Actions.
- Arbitrary-precision checks use Arb through `python-flint`; the baseline is **160 bits**, with higher-precision reruns to detect numerical artifacts.
- Computation is evidence and falsification support, not a proof of regularity.

## Current physical hypothesis

For a coherent structure with characteristic transverse area `A`, length `L`, circulation `Gamma`, and viscosity `nu`, the first reduced model gives

- viscous power scale: `P_nu ~ nu * Gamma^2 * L / A`,
- nonlinear turnover time: `tau_nl ~ A / Gamma`,
- one-turnover viscous toll: `P_nu * tau_nl ~ nu * Gamma * L`.

The striking feature is that the transverse area cancels. The first computational program asks whether this scale-independent cancellation is robust under controlled deformations of the model, and then progressively replaces modeling assumptions by quantities extracted directly from the PDE.

## Experiments

`00_precision_and_scaling.py` is only a numerical/symbolic sanity gate. It checks that the proposed scale cancellation is not a floating-point artifact and records its stability at several Arb precisions.
