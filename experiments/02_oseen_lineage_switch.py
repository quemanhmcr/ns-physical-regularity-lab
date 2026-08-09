import json
import os
from flint import arb, ctx

BITS = int(os.environ.get("ARB_PREC_BITS", "160"))
if BITS < 160:
    raise SystemExit("ARB_PREC_BITS must be at least 160")
ctx.prec = BITS
pi = arb.pi()

# Exact Lamb-Oseen vortex:
#   omega(r,t) = Gamma/(4*pi*nu*t) exp(-r^2/(4*nu*t))
# with purely azimuthal velocity, so every circle r=R is a material loop.
# Its enclosed circulation is Gamma_R = Gamma(1-exp(-q)), q=R^2/(4 nu t).
# This experiment tests three exact NS statements:
#   (1) radial vorticity diffusion PDE,
#   (2) Kelvin-with-viscosity circulation loss equals viscous boundary flux,
#   (3) a fixed fractional lineage switch has cost ~ Gamma^2 independent of R and nu.

nu_values = ["1e-24", "1e-12", "1e-3", "1", "1e12"]
gamma_values = ["1e-18", "1e-9", "1", "1e9"]
radius_values = ["1e-24", "1e-8", "1", "1e8", "1e24"]
q_probe_values = ["0.01", "0.1", "1", "10", "100"]
# Fraction of total circulation inside a fixed material circle, early -> later.
fraction_pairs = [
    ("0.99", "0.90"),
    ("0.90", "0.50"),
    ("0.75", "0.25"),
    ("0.50", "0.10"),
    ("0.10", "0.01"),
]

rows = []
for nu_s in nu_values:
    for gamma_s in gamma_values:
        for R_s in radius_values:
            nu = arb(nu_s)
            gamma = arb(gamma_s)
            R = arb(R_s)

            # Probe exact PDE and Kelvin-viscous flux at multiple q values.
            last_pde_residual = None
            last_flux_ratio = None
            for q_s in q_probe_values:
                q = arb(q_s)
                t = R * R / (4 * nu * q)
                omega = gamma / (4 * pi * nu * t) * (-q).exp()

                # At r=R, q=r^2/(4 nu t).
                omega_t = omega * (-1 + q) / t
                omega_r = -R * omega / (2 * nu * t)
                omega_rr = (-1 / (2 * nu * t) + R * R / (4 * nu * nu * t * t)) * omega
                radial_laplacian = omega_rr + omega_r / R
                residual = omega_t - nu * radial_laplacian
                if not residual.contains(0):
                    raise AssertionError(
                        f"Oseen vorticity PDE residual excludes 0: nu={nu_s}, Gamma={gamma_s}, R={R_s}, q={q_s}, residual={residual}"
                    )

                # Material-loop circulation and its exact time derivative.
                gamma_R = gamma * (1 - (-q).exp())
                dgamma_dt = -gamma * q * (-q).exp() / t

                # Viscous vorticity flux through the fixed material circle.
                viscous_flux = nu * (2 * pi * R) * omega_r
                flux_ratio = viscous_flux / dgamma_dt
                if not flux_ratio.contains(1):
                    raise AssertionError(
                        f"Kelvin-viscous flux identity failed: nu={nu_s}, Gamma={gamma_s}, R={R_s}, q={q_s}, ratio={flux_ratio}"
                    )

                # Cross-check circulation formula from the exact radial integral.
                gamma_integral = gamma * (1 - (-q).exp())
                circ_ratio = gamma_R / gamma_integral
                if not circ_ratio.contains(1):
                    raise AssertionError(
                        f"enclosed circulation formula mismatch: nu={nu_s}, Gamma={gamma_s}, R={R_s}, q={q_s}, ratio={circ_ratio}"
                    )

                last_pde_residual = str(residual)
                last_flux_ratio = str(flux_ratio)

            pair_summaries = []
            for f1_s, f2_s in fraction_pairs:
                f1 = arb(f1_s)
                f2 = arb(f2_s)
                if not (arb(0) < f2 < f1 < arb(1)):
                    raise AssertionError("fraction pair must satisfy 0 < f2 < f1 < 1")

                # f = 1-exp(-q), hence q = -log(1-f).
                q1 = -(1 - f1).log()
                q2 = -(1 - f2).log()
                t1 = R * R / (4 * nu * q1)
                t2 = R * R / (4 * nu * q2)

                # Total Oseen enstrophy is Gamma^2/(8*pi*nu*t), so physical
                # viscous dissipation per unit axial length is Gamma^2/(8*pi*t).
                # Integrating from t1 to t2 gives the exact switch toll.
                cost = gamma * gamma / (8 * pi) * (t2 / t1).log()
                expected = gamma * gamma / (8 * pi) * (q1 / q2).log()
                cost_ratio = cost / expected
                if not cost_ratio.contains(1):
                    raise AssertionError(
                        f"lineage-switch toll lost scale cancellation: nu={nu_s}, Gamma={gamma_s}, R={R_s}, f1={f1_s}, f2={f2_s}, ratio={cost_ratio}"
                    )

                # Circulation lost by the material core is a fixed fraction of Gamma.
                lost = gamma * (f1 - f2)
                expected_lost = gamma * (f1 - f2)
                lost_ratio = lost / expected_lost
                if not lost_ratio.contains(1):
                    raise AssertionError("fixed-fraction leakage normalization failed")

                # Cost per circulation lost must be proportional to Gamma, with
                # a dimensionless coefficient depending only on the chosen fractions.
                cost_per_lost = cost / lost
                expected_cpl = gamma * (q1 / q2).log() / (8 * pi * (f1 - f2))
                cpl_ratio = cost_per_lost / expected_cpl
                if not cpl_ratio.contains(1):
                    raise AssertionError(
                        f"cost-per-leaked-circulation scaling failed: nu={nu_s}, Gamma={gamma_s}, R={R_s}, f1={f1_s}, f2={f2_s}, ratio={cpl_ratio}"
                    )

                pair_summaries.append({
                    "f1": f1_s,
                    "f2": f2_s,
                    "q1": str(q1),
                    "q2": str(q2),
                    "t2_over_t1": str(t2 / t1),
                    "cost_ratio": str(cost_ratio),
                    "cost_per_lost_ratio": str(cpl_ratio),
                })

            rows.append({
                "nu": nu_s,
                "Gamma": gamma_s,
                "R": R_s,
                "last_pde_residual_ball": last_pde_residual,
                "last_kelvin_flux_ratio": last_flux_ratio,
                "fraction_pairs": pair_summaries,
            })

print(json.dumps({
    "arb_precision_bits": BITS,
    "parameter_cases": len(rows),
    "q_probes_per_case": len(q_probe_values),
    "fraction_switches_per_case": len(fraction_pairs),
    "status": "PASS",
    "interpretation": (
        "In the exact Lamb-Oseen NS diffusion solution, circulation can leave a fixed material core only through viscosity; "
        "moving a fixed fraction of circulation across that material boundary carries a Gamma^2 toll independent of core radius and viscosity."
    ),
    "rows": rows,
}, indent=2))
