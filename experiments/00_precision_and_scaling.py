import json
import os
from flint import arb, ctx
import mpmath as mp
from fractions import Fraction

BITS = int(os.environ.get("ARB_PREC_BITS", "160"))
if BITS < 160:
    raise SystemExit("ARB_PREC_BITS must be at least 160")
ctx.prec = BITS
mp.mp.dps = max(80, int(BITS * 0.30103) + 30)

nu_s = "0.000001"
gamma_s = "0.000009"
L_s = "3.125"

nu = arb(nu_s)
gamma = arb(gamma_s)
L = arb(L_s)

nu_mp = mp.mpf(nu_s)
gamma_mp = mp.mpf(gamma_s)
L_mp = mp.mpf(L_s)
expected_mp = nu_mp * gamma_mp * L_mp

rows = []
for exponent in range(-60, 61, 10):
    A_s = f"1e{exponent}"
    A = arb(A_s)
    P = nu * gamma * gamma * L / A
    tau_nl = A / gamma
    toll = P * tau_nl

    A_mp = mp.mpf(A_s)
    P_mp = nu_mp * gamma_mp * gamma_mp * L_mp / A_mp
    tau_mp = A_mp / gamma_mp
    toll_mp = P_mp * tau_mp
    rel_mp = abs((toll_mp - expected_mp) / expected_mp)

    rows.append({
        "A": A_s,
        "arb_toll": str(toll),
        "mp_toll": mp.nstr(toll_mp, 70),
        "mp_relative_error": mp.nstr(rel_mp, 20),
    })

    if rel_mp != 0:
        raise AssertionError(f"mpmath cancellation failed at A={A_s}: {rel_mp}")

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
    "expected_toll_mpmath": mp.nstr(expected_mp, 70),
    "scales_tested": len(rows),
    "status": "PASS",
    "rows": rows,
}
print(json.dumps(summary, indent=2))
