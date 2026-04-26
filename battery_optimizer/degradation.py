"""
degradation.py
--------------
Degradation cost module for the BESS aging-aware EMS optimization.

Two formulations are provided:

1) LINEAR (throughput-based):
       C_deg = v_SoH * alpha * (c + d) * dt
   Matches the "objective function first draft" of the report.
   Easy to plug into an LP/MILP, but does NOT capture the DoD non-linearity.

2) WÖHLER / MINER'S RULE (DoD-based, recommended):
       N(DoD) = a * DoD^(-b)            (Wöhler curve)
       wear_per_cycle = 1 / N(DoD)      (Miner's rule)
       C_deg = C_replacement * wear_per_cycle
   Captures the key insight from the literature: deep cycles are
   disproportionately damaging compared to shallow ones.

You can call either one, or use `degradation_cost(..., model="wohler")`
as the single entry point that the rest of the codebase uses.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence, Literal, Union
import numpy as np
import pandas as pd


# ----------------------------------------------------------------------------
# Default location of the data files (same directory as this module)
# ----------------------------------------------------------------------------
_THIS_DIR = Path(__file__).resolve().parent
DEFAULT_PARAMS_CSV = _THIS_DIR / "battery_params.csv"
DEFAULT_CYCLES_CSV = _THIS_DIR / "cycle_life_data.csv"


# ----------------------------------------------------------------------------
# Data class for battery parameters
# ----------------------------------------------------------------------------
@dataclass
class BatteryParams:
    """Per-container battery parameters needed for degradation costing."""
    E_max: float                              # usable energy capacity (MWh)
    eta_c: float                              # charging efficiency
    eta_d: float                              # discharging efficiency
    a: float                                  # Wöhler scale parameter
    b: float                                  # Wöhler shape parameter
    replacement_cost_per_MWh: float           # $ / MWh of capacity
    # Linear-model coefficients (auto-derived from the Wöhler params)
    alpha_lin: float = field(default=None)
    v_SoH: float = field(default=None)

    def __post_init__(self):
        # Auto-derive linear model coefficients so the two models are
        # roughly comparable at the reference DoD (80% by convention).
        if self.alpha_lin is None:
            N_eq = self.a * (0.8) ** (-self.b)
            self.alpha_lin = 1.0 / (2.0 * N_eq * self.E_max)
        if self.v_SoH is None:
            # 1 unit of SoH lost = full replacement of the energy capacity
            self.v_SoH = self.replacement_cost_per_MWh * self.E_max

    # ------------------------------------------------------------------
    # CSV loaders
    # ------------------------------------------------------------------
    @classmethod
    def from_csv(cls,
                 params_csv: Union[str, Path] = DEFAULT_PARAMS_CSV,
                 cycles_csv: Union[str, Path] = None,
                 refit_wohler: bool = False) -> "BatteryParams":
        """
        Load battery parameters from `battery_params.csv`.

        If `refit_wohler=True`, the Wöhler coefficients (a, b) are recomputed
        from `cycle_life_data.csv` instead of being read from the params file.
        Useful if you edit the cycle-life CSV and want the params to follow.

        Parameters
        ----------
        params_csv : path to the scalar-parameters CSV.
        cycles_csv : path to the cycle-life data CSV (only used when refitting).
        refit_wohler : whether to re-run the regression on cycle_life_data.csv.
        """
        df = pd.read_csv(params_csv)
        # turn the (parameter, value) long format into a dict
        p = dict(zip(df["parameter"], df["value"].astype(float)))

        if refit_wohler:
            cycles_csv = cycles_csv or DEFAULT_CYCLES_CSV
            a, b = fit_wohler_from_csv(cycles_csv)
        else:
            a, b = p["wohler_a"], p["wohler_b"]

        return cls(
            E_max=p["E_max"],
            eta_c=p["eta_c"],
            eta_d=p["eta_d"],
            a=a,
            b=b,
            replacement_cost_per_MWh=p["replacement_cost_per_MWh"],
        )


# ----------------------------------------------------------------------------
# Helper: refit the Wöhler curve from the cycle-life CSV
# ----------------------------------------------------------------------------
def fit_wohler_from_csv(cycles_csv: Union[str, Path] = DEFAULT_CYCLES_CSV
                        ) -> tuple[float, float]:
    """
    Fit N(DoD) = a * DoD^(-b) on the (DoD, cycles) points stored in the
    cycle-life CSV and return (a, b).
    """
    from scipy.optimize import curve_fit
    df = pd.read_csv(cycles_csv)
    DoD = df["DoD_fraction"].to_numpy(dtype=float)
    N = df["cycles_to_EOL"].to_numpy(dtype=float)
    (a, b), _ = curve_fit(lambda d, a, b: a * d ** (-b), DoD, N, p0=[4000.0, 1.0])
    return float(a), float(b)


# ----------------------------------------------------------------------------
# Helper: extract half-cycles via rainflow-lite (full discharge events)
# ----------------------------------------------------------------------------
def _detect_discharge_events(soc: np.ndarray):
    """
    Walk the SoC trajectory and return the list of (DoD) of each
    discharge half-cycle. A discharge half-cycle is a monotonic
    decrease in SoC between a local max and the next local min.

    This is a simplified rainflow that's adequate for daily BESS
    schedules where cycles are usually well-defined.
    """
    soc = np.asarray(soc, dtype=float)
    if len(soc) < 2:
        return []

    # Find turning points (local extrema, including endpoints)
    diffs = np.diff(soc)
    sign_changes = np.diff(np.sign(diffs)) != 0
    turning_idx = np.concatenate(([0], np.where(sign_changes)[0] + 1, [len(soc) - 1]))
    turning_soc = soc[turning_idx]

    # DoD = drop from a local max to the next local min
    DoDs = []
    for i in range(len(turning_soc) - 1):
        delta = turning_soc[i] - turning_soc[i + 1]
        if delta > 0:                   # discharge event
            DoDs.append(delta)
    return DoDs


# ----------------------------------------------------------------------------
# Model 1 -- Linear (throughput proportional)
# ----------------------------------------------------------------------------
def degradation_cost_linear(c: Sequence[float],
                            d: Sequence[float],
                            params: BatteryParams,
                            dt: float = 1.0) -> float:
    """
    C_deg = v_SoH * alpha * sum_t (c_t + d_t) * dt

    Parameters
    ----------
    c, d : sequences of charging / discharging power (MW), length T
    params : BatteryParams
    dt : timestep duration (h)
    """
    c = np.asarray(c, dtype=float)
    d = np.asarray(d, dtype=float)
    throughput = np.sum(c + d) * dt                       # MWh moved
    return params.v_SoH * params.alpha_lin * throughput   # $


# ----------------------------------------------------------------------------
# Model 2 -- Wöhler / Miner's rule (DoD aware)  -- RECOMMENDED
# ----------------------------------------------------------------------------
def degradation_cost_wohler(c: Sequence[float],
                            d: Sequence[float],
                            params: BatteryParams,
                            dt: float = 1.0,
                            e0: float = None) -> float:
    """
    Compute degradation cost using the Wöhler/Miner formulation.

    Steps:
      1. Reconstruct the SoC trajectory from (c, d) and initial energy e0.
      2. Detect discharge half-cycles and their DoD.
      3. For each DoD: wear = 0.5 / N(DoD)   (half-cycle, hence the 0.5)
                       cost = wear * total_replacement_cost
      4. Sum across all events.

    Parameters
    ----------
    c, d : charging / discharging power (MW)
    params : BatteryParams
    dt : timestep (h)
    e0 : initial stored energy (MWh). Defaults to 50% SoC.
    """
    c = np.asarray(c, dtype=float)
    d = np.asarray(d, dtype=float)
    T = len(c)
    if e0 is None:
        e0 = 0.5 * params.E_max

    # Energy trajectory in MWh, then convert to SoC fraction
    e = np.empty(T + 1)
    e[0] = e0
    for t in range(T):
        e[t + 1] = e[t] + params.eta_c * c[t] * dt - (1.0 / params.eta_d) * d[t] * dt
    soc = e / params.E_max

    DoDs = _detect_discharge_events(soc)

    # Total $ if the battery were fully replaced (1 unit SoH lost)
    total_replacement_cost = params.replacement_cost_per_MWh * params.E_max

    cost = 0.0
    for DoD in DoDs:
        DoD = max(DoD, 1e-3)                              # numerical floor
        N = params.a * DoD ** (-params.b)                 # cycles to EOL
        wear = 0.5 / N                                    # half-cycle (Miner)
        cost += wear * total_replacement_cost
    return cost


# ----------------------------------------------------------------------------
# Unified entry point used by total_cost(...)
# ----------------------------------------------------------------------------
def degradation_cost(p_k,
                     battery_params: BatteryParams,
                     dt: float = 1.0,
                     model: Literal["linear", "wohler"] = "wohler",
                     e0: float = None) -> float:
    """
    Drop-in replacement for the stub in your skeleton.

    `p_k` may be either:
      - a 1-D array of net power injection per container (positive = discharge),
        in which case c = max(-p,0), d = max(p,0).
      - a tuple (c, d) of two arrays already separated.
    """
    if isinstance(p_k, tuple) and len(p_k) == 2:
        c, d = p_k
    else:
        p = np.asarray(p_k, dtype=float)
        c = np.maximum(-p, 0.0)
        d = np.maximum(p, 0.0)

    if model == "linear":
        return degradation_cost_linear(c, d, battery_params, dt)
    elif model == "wohler":
        return degradation_cost_wohler(c, d, battery_params, dt, e0)
    else:
        raise ValueError(f"Unknown model: {model!r}")


# ----------------------------------------------------------------------------
# Quick self-test / sanity check
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    # Load all parameters from CSV (no hardcoded values)
    bp = BatteryParams.from_csv()

    print("=" * 60)
    print("Loaded from battery_params.csv:")
    print(f"  E_max = {bp.E_max} MWh")
    print(f"  eta_c = {bp.eta_c}, eta_d = {bp.eta_d}")
    print(f"  Wöhler a = {bp.a:.2f}, b = {bp.b:.4f}")
    print(f"  Replacement = ${bp.replacement_cost_per_MWh:,.0f} / MWh")
    print(f"  v_SoH = ${bp.v_SoH:,.0f}, alpha_lin = {bp.alpha_lin:.2e}")
    print()
    print("=" * 60)
    print("Calibration check (Wöhler curve from CSV)")
    print("=" * 60)
    for DoD_pct in [50, 70, 80, 100]:
        DoD = DoD_pct / 100
        N = bp.a * DoD ** (-bp.b)
        print(f"DoD = {DoD_pct:3d}%  ->  N = {N:6.0f} cycles to EOL")

    print()
    print("=" * 60)
    print("Scenario A: one full daily cycle (DoD = 70%, operational range)")
    print("=" * 60)
    # 12h charge at 0.5 * P_max, 12h idle, then discharge --> simplified
    T = 24
    c = np.zeros(T); d = np.zeros(T)
    # Charge 15% -> 85% over 4h (operational range), then discharge over 4h
    # Capacity is 3.916 MWh, 70% of that = 2.74 MWh, over 4h = 0.685 MW
    c[2:6] = 0.685   # charge 02:00->06:00
    d[18:22] = 0.685 # discharge 18:00->22:00, e0=15% SoC

    e0 = 0.15 * bp.E_max
    print(f"  Linear  : ${degradation_cost((c, d), bp, model='linear'):,.2f}")
    print(f"  Wöhler  : ${degradation_cost((c, d), bp, model='wohler', e0=e0):,.2f}")

    print()
    print("=" * 60)
    print("Scenario B: same throughput, but split into 2 shallow cycles")
    print("=" * 60)
    c = np.zeros(T); d = np.zeros(T)
    # Two shallow cycles: half the energy each
    c[2:4] = 0.685; d[8:10] = 0.685       # cycle 1 (DoD ~35%)
    c[12:14] = 0.685; d[18:20] = 0.685    # cycle 2 (DoD ~35%)
    e0 = 0.50 * bp.E_max
    print(f"  Linear  : ${degradation_cost((c, d), bp, model='linear'):,.2f}")
    print(f"  Wöhler  : ${degradation_cost((c, d), bp, model='wohler', e0=e0):,.2f}")
    print()
    print("Notice: linear sees the SAME cost (same throughput).")
    print("        Wöhler rewards the 2-shallow-cycle strategy: this is the")
    print("        non-linearity the literature emphasizes.")
