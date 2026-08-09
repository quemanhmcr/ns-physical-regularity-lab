import json
import os
from fractions import Fraction
from flint import arb, ctx
import mpmath as mp

BITS = int(os.environ.get("ARB_PREC_BITS", "160"))
if BITS < 160:
    raise SystemExit("ARB_PREC_BITS must be at least 160")
ctx.prec = BITS
mp.mp.dps = max(100, int(BITS * 0.30103) + 40)

nu_s = "0.000001"
gamma_s = "0.000009"
L_s = "3.125"

nu = arb(nu_s)
gamma = arb(gamma_s)
L = arb(L_s)
expected_arb = nu * gamma * L

nu_mp = mp.mpf(nu_s)
gamma_mp = mp.mpf(gamma_s)
L_mp = mp.mpf(L_s)
expected_mp = nu_mp * gamma_mp * L_mp
# Demand roughly BITS correct binary digits in the cancellation test, while
# retaining guard digits in the mpmath working precision.
mp_rel_tol = mp.mpf(2) ** (-BITS)

rows = []
for exponent in range(-60, 61, 10):
    A_s = f"1e{exponent}"
    A = arb(A_s)
    P = nu * gamma * gamma * L / A
    tau_nl = A / gamma
    toll = P * tau_nl
    ratio_arb = toll / expected_arb

    # Arb is an interval/ball computation.  The mathematically exact value 1
    # must be enclosed by the computed ratio at every tested scale.
    if not ratio_arb.contains(1):
        raise AssertionError(f"Arb enclosure lost exact ratio 1 at A={A_s}: {ratio_arb}")

    A_mp = mp.mpf(A_s)
    P_mp = nu_mp * gamma_mp * gamma_mp * L_mp / A_mp
    tau_mp = A_mp / gamma_mp
    toll_mp = P_mp * tau_mp
    rel_mp = abs((toll_mp - expected_mp) / expected_mp)
    if rel_mp > mp_rel_tol:
        raise AssertionError(
            f"mpmath relative error too large at A={A_s}: {rel_mp} > {mp_rel_tol}"
        )

    rows.append({
        "A": A_s,
        "arb_ratio": str(ratio_arb),
        "arb_toll": str(toll),
        "mp_toll": mp.nstr(toll_mp, 80),
        "mp_relative_error": mp.nstr(rel_mp, 30),
    })

# Independent exact arithmetic gate.  This proves only the algebraic reduced
# identity, not its physical assumptions.
nu_q = Fraction(1, 10**6)
gamma_q = Fraction(9, 10**6)
L_q = Fraction(25, 8)
for exponent in range(-30, 31, 5):
    A_q = Fraction(10**exponent, 1) if exponent >= 0 else Fraction(1, 10**(-exponent))
    P_q = nu_q * gamma_q * gamma_q * L_q / A_q
    tau_q = A_q / gamma_q
    ratio_q = (P_q * tau_q) / (nu_q * gamma_q * L_q)
    if ratio_q != 1:
        raise AssertionError(f"exact rational cancellation failed at exponent {exponent}")

summary = {
    "arb_precision_bits": BITS,
    "mp_relative_tolerance": mp.nstr(mp_rel_tol, 30),
    "expected_toll_mpmath": mp.nstr(expected_mp, 80),
    "scales_tested": len(rows),
    "status": "PASS",
    "rows": rows,
}
print(json.dumps(summary, indent=2))
