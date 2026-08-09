import json
import os
from flint import arb, ctx

BITS = int(os.environ.get("ARB_PREC_BITS", "160"))
if BITS < 160:
    raise SystemExit("ARB_PREC_BITS must be at least 160")
ctx.prec = BITS
pi = arb.pi()

# Sweep several decades.  All checks are interval checks with Arb: the exact
# target value must lie inside the computed ball.
a_values = ["1e-24", "1e-12", "1", "1e12", "1e24"]
nu_values = ["1e-18", "1e-9", "1e-3", "1"]
gamma_values = ["1e-12", "1e-6", "1", "1e6"]
radius_multipliers = ["0.125", "0.5", "1", "2", "8"]

rows = []
for a_s in a_values:
    for nu_s in nu_values:
        for gamma_s in gamma_values:
            a = arb(a_s)
            nu = arb(nu_s)
            gamma = arb(gamma_s)

            beta = a / (4 * nu)
            delta2 = 4 * nu / a
            delta = delta2.sqrt()
            omega0 = gamma * a / (4 * pi * nu)

            # Exact Gaussian integrals for this profile.
            circulation = omega0 * (4 * pi * nu / a)
            circulation_ratio = circulation / gamma
            if not circulation_ratio.contains(1):
                raise AssertionError(
                    f"circulation normalization failed: a={a_s}, nu={nu_s}, Gamma={gamma_s}, ratio={circulation_ratio}"
                )

            enstrophy_per_length = gamma * gamma * a / (8 * pi * nu)
            viscous_vorticity_power_per_length = nu * enstrophy_per_length
            strain_time = 1 / a
            toll = viscous_vorticity_power_per_length * strain_time
            expected_toll = gamma * gamma / (8 * pi)
            toll_ratio = toll / expected_toll
            if not toll_ratio.contains(1):
                raise AssertionError(
                    f"scale-independent toll failed: a={a_s}, nu={nu_s}, Gamma={gamma_s}, ratio={toll_ratio}"
                )

            # Verify the steady axisymmetric vorticity equation pointwise:
            #   -(a r/2) omega_r = a omega + nu(omega_rr + omega_r/r)
            # for the Gaussian profile.  Derivatives are inserted analytically;
            # Arb then checks that the residual ball encloses zero over a wide
            # range of radii relative to the dynamically generated core radius.
            worst_residual = None
            for m_s in radius_multipliers:
                r = arb(m_s) * delta
                omega = omega0 * (-beta * r * r).exp()
                omega_r = -2 * beta * r * omega
                omega_rr = (-2 * beta + 4 * beta * beta * r * r) * omega
                lap_radial = omega_rr + omega_r / r
                residual = (-a * r * omega_r / 2) - (a * omega + nu * lap_radial)
                if not residual.contains(0):
                    raise AssertionError(
                        f"NS vorticity residual excludes 0: a={a_s}, nu={nu_s}, Gamma={gamma_s}, r/delta={m_s}, residual={residual}"
                    )
                worst_residual = str(residual)

            rows.append({
                "a": a_s,
                "nu": nu_s,
                "Gamma": gamma_s,
                "delta2": str(delta2),
                "circulation_ratio": str(circulation_ratio),
                "toll_ratio": str(toll_ratio),
                "last_residual_ball": worst_residual,
            })

print(json.dumps({
    "arb_precision_bits": BITS,
    "cases": len(rows),
    "radius_checks_per_case": len(radius_multipliers),
    "status": "PASS",
    "interpretation": "Exact strained-vortex NS balance preserves a circulation-squared toll per strain time independent of core scale.",
    "rows": rows,
}, indent=2))
