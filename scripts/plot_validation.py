#!/usr/bin/env python3
"""
Three-panel validation plot for BBRsim crack dispatch.

Plots (from build directory):
  conda run -n bbrsim python ../scripts/plot_validation.py crack1_output.csv crack2_output.csv

Panel 1 — Exit z-position: crack1 clusters at z≈0, crack2 at z≈3 mm.
Panel 2 — Transmitted k_z: crack2 (wider gap) → narrower angular spread.
Panel 3 — Transmittance: simulated T vs HFSS table value per crack.
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ---------------------------------------------------------------------------

def hfss_T(dataset_id, hfss_root, iwaveTheta=180., iwavePhi=0.):
    """Return clamped expected transmittance for equal-mix E_theta/E_phi polarization."""
    T_by_ephi = {}
    for ephi in (0, 1):
        csv = os.path.join(hfss_root, f"{dataset_id}_Ephi={ephi}", "waveguide.csv")
        df  = pd.read_csv(csv)
        row = df[
            (np.abs(df["IWaveTheta"] - iwaveTheta) < 1e-5) &
            (np.abs(df["IWavePhi"]   - iwavePhi)   < 1e-5)
        ].iloc[0]
        T_by_ephi[ephi] = row["OutgoingPower"] / row["IngoingPower"]
    # gun fires equal E_theta / E_phi → 50/50 weight
    return min(0.5 * T_by_ephi[0] + 0.5 * T_by_ephi[1], 1.0)

# ---------------------------------------------------------------------------

def main():
    f1 = sys.argv[1] if len(sys.argv) > 1 else "crack1_output.csv"
    f2 = sys.argv[2] if len(sys.argv) > 2 else "crack2_output.csv"

    df1 = pd.read_csv(f1)
    df2 = pd.read_csv(f2)

    id1 = df1["crack_id"].iloc[0]
    id2 = df2["crack_id"].iloc[0]

    script_dir = os.path.dirname(os.path.abspath(__file__))
    hfss_root  = os.path.join(script_dir, "..", "HFSSSimData")

    T_exp1 = hfss_T(id1, hfss_root)
    T_exp2 = hfss_T(id2, hfss_root)

    t1 = df1[df1["transmitted"] == 1]
    t2 = df2[df2["transmitted"] == 1]

    C1, C2 = "steelblue", "tomato"

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        "BBRsim crack dispatch validation — 500 GHz, normal incidence",
        fontsize=13, y=1.01
    )

    # ------------------------------------------------------------------
    # Panel 1: exit z-position (transmitted only)
    # ------------------------------------------------------------------
    ax = axes[0]
    z1_mm = t1["pos_z_m"].dropna() * 1e3
    z2_mm = t2["pos_z_m"].dropna() * 1e3
    lo = min(z1_mm.min() if len(z1_mm) else -1., z2_mm.min() if len(z2_mm) else -1.) - 1.
    hi = max(z1_mm.max() if len(z1_mm) else 4.,  z2_mm.max() if len(z2_mm) else 4.)  + 1.
    bins = np.linspace(lo, hi, 60)
    ax.hist(z1_mm, bins=bins, alpha=0.65, color=C1, label="crack1  b=50 µm  (z_gun=0 mm)")
    ax.hist(z2_mm, bins=bins, alpha=0.65, color=C2, label="crack2  b=100 µm  (z_gun=3 mm)")
    ax.axvline(0., color=C1, linestyle="--", linewidth=1.4)
    ax.axvline(3., color=C2, linestyle="--", linewidth=1.4)
    ax.set_xlabel("Exit z-position (mm)")
    ax.set_ylabel("Transmitted events")
    ax.set_title("Panel 1 — Exit z-position\nEach crack exits near its own center")
    ax.legend(fontsize=8)

    # ------------------------------------------------------------------
    # Panel 2: transmitted k_z distribution
    # ------------------------------------------------------------------
    ax = axes[1]
    ax.hist(t1["dir_z"].dropna(), bins=60, alpha=0.65, color=C1,
            density=True, label="crack1  b=50 µm")
    ax.hist(t2["dir_z"].dropna(), bins=60, alpha=0.65, color=C2,
            density=True, label="crack2  b=100 µm")
    ax.set_xlabel(r"Transmitted $k_z$  (unit vector component)")
    ax.set_ylabel("Density")
    ax.set_title("Panel 2 — Diffracted $k_z$ spread\nWider gap → narrower angular spread")
    ax.legend(fontsize=8)

    # ------------------------------------------------------------------
    # Panel 3: transmittance — simulated vs HFSS
    # ------------------------------------------------------------------
    ax = axes[2]
    T_sim1 = (df1["transmitted"] == 1).mean()
    T_sim2 = (df2["transmitted"] == 1).mean()

    x = np.array([0., 1.])
    w = 0.3
    b_sim  = ax.bar(x - w/2, [T_sim1, T_sim2], w,
                    color=[C1, C2], alpha=0.85, label="Simulated")
    b_exp  = ax.bar(x + w/2, [T_exp1, T_exp2], w,
                    color=[C1, C2], alpha=0.4, hatch="//", label="HFSS expected")

    for bar, v in zip(b_sim, [T_sim1, T_sim2]):
        ax.text(bar.get_x() + bar.get_width() / 2., v + 0.012,
                f"{v:.1%}", ha="center", va="bottom", fontsize=9, fontweight="bold")
    for i, v in enumerate([T_exp1, T_exp2]):
        ax.text(x[i] + w/2., v + 0.012,
                f"{v:.1%}", ha="center", va="bottom", fontsize=9, color="gray")

    n1, n2 = len(df1), len(df2)
    ax.set_xticks(x)
    ax.set_xticklabels([f"crack1\nb=50 µm\n(N={n1})", f"crack2\nb=100 µm\n(N={n2})"])
    ax.set_ylabel("Transmittance")
    ax.set_ylim(0., 1.)
    ax.set_title("Panel 3 — Transmittance\nSimulated vs HFSS table")
    ax.legend(fontsize=8)

    # ------------------------------------------------------------------
    plt.tight_layout()
    out = "validation_plots.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")
    plt.show()


if __name__ == "__main__":
    main()
