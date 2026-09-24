"""
plot_crack_frequency.py
Visualise the frequency-keyed HFSS lookup on the output of crack_frequency.mac
(16 runs against the mock five-frequency tree). Four panels:

  a) Selection map   - which grid dataset each photon was given, against its own
                       frequency. The staircase comes from the broadband Planck
                       run; the fixed-gun probes are marked on top. Dashed lines
                       are the log-space midpoints where the choice switches.
  b) Transmittance   - observed per gun run against the prediction computed from
                       the mock dataset that run selected. Agreement here is what
                       proves the recorded frequency is the dataset actually
                       sampled, since the mock scales T per frequency.
  c) Planck spectrum - the 20 K photon-number spectrum entering the crack,
                       coloured by the dataset each photon received.
  d) Exit signature  - fraction of transmitted photons leaving with k_z < 0, for
                       the normal-incidence gun runs. The 150 and 1500 GHz mock
                       sets have their far field truncated to Theta <= 90 deg, so
                       they read zero; the others sit near one half. Only normal
                       incidence shows this: at an oblique azimuth the
                       quarter-symmetry unfold mirrors phi_hat for half the
                       photons, which flips the sign of k_z in the world frame.

Colour: the five grid frequencies are ordered levels, so they use one blue
ordinal ramp (light = low frequency), not categorical hues. Steps 250-700 of
the reference ramp; the lightest clears the surface at 2.06:1.

Usage:
    conda run -n bbrsim python scripts/plot_crack_frequency.py [output_dir] \\
        [--data-dir MOCK_ROOT] [--out PATH]
    defaults: build/output, <output_dir>/../../mock_hfss,
              <output_dir>/crack_frequency_overview.png
"""
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim import hfss, select
from bbrsim.io import load_many

GRID = [50.0, 150.0, 500.0, 1500.0, 5000.0]
MIDS = [np.sqrt(a * b) for a, b in zip(GRID[:-1], GRID[1:])]
CRACK1 = "InfParallelPlate_crack1Rohan"
TRUNCATED = {150.0, 1500.0}          # mock sets with far field cut to Theta <= 90 deg
PLANCK_RUN = 15

# One blue ordinal ramp, light -> dark with frequency (validated: monotone
# lightness, adjacent dL >= 0.06, light end 2.06:1 on the surface, single hue).
RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]
COLOR = dict(zip(GRID, RAMP))
INK, INK_2, INK_MUTED = "#0b0b0b", "#52514e", "#8a8a85"
SURFACE = "#fcfcfb"

ap = argparse.ArgumentParser()
ap.add_argument("output_dir", nargs="?", default="build/output")
ap.add_argument("--data-dir", default=None, help="mock data ROOT (contains waveguides/)")
ap.add_argument("--out", default=None)
args = ap.parse_args()
mock_root = args.data_dir or os.path.normpath(
    os.path.join(args.output_dir, "..", "..", "mock_hfss"))
wg_dir = os.path.join(mock_root, "waveguides")
out_path = args.out or os.path.join(args.output_dir, "crack_frequency_overview.png")

files = sorted(glob.glob(os.path.join(args.output_dir, "bbr_freq_r*.root")))
if not files:
    sys.exit(f"no bbr_freq_r*.root under {args.output_dir}")

df = load_many(files)
entries = select.crack_crossings(df)
decided = entries[entries["status"].isin(
    ["BBRDiffractionTransmit", "BBRDiffractionReflect"])]
crack1 = decided[decided["vol_post"].str.contains("crack1", na=False)].copy()
crack1["nu_GHz"] = hfss.photon_frequency_GHz(crack1["energy_eV"].to_numpy(float))

planck = crack1[crack1["run_id"] == PLANCK_RUN]
gun = crack1[crack1["run_id"] != PLANCK_RUN]

fig, axes = plt.subplots(2, 2, figsize=(14, 9), facecolor=SURFACE)
for ax in axes.flat:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_MUTED)
    ax.tick_params(colors=INK_2, labelsize=9)
    ax.grid(True, which="major", color=INK_MUTED, alpha=0.18, linewidth=0.6)
    ax.set_axisbelow(True)

# --- (a) selection map -------------------------------------------------------
ax = axes[0, 0]
for f in GRID:
    sub = planck[planck["hfss_freq_GHz"] == f]
    ax.scatter(sub["nu_GHz"], sub["hfss_freq_GHz"], s=5, color=COLOR[f],
               alpha=0.5, linewidths=0, zorder=2)
for m in MIDS:
    ax.axvline(m, color=INK_MUTED, linestyle="--", linewidth=0.9, zorder=1)
gun_pts = gun.groupby("run_id").agg(nu=("nu_GHz", "first"),
                                    sel=("hfss_freq_GHz", "first"))
ax.scatter(gun_pts["nu"], gun_pts["sel"], s=70, facecolors="none",
           edgecolors=INK, linewidths=1.4, zorder=3)
for label, x in (("clamped\nlow", 20.0), ("clamped\nhigh", 10000.0)):
    row = gun_pts[np.isclose(gun_pts["nu"], x, rtol=1e-3)]
    if len(row):
        ax.annotate(label, (row["nu"].iloc[0], row["sel"].iloc[0]),
                    textcoords="offset points", xytext=(0, -34), ha="center",
                    fontsize=8, color=INK_2)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_yticks(GRID); ax.set_yticklabels([f"{f:g}" for f in GRID])
ax.set_xlabel("photon frequency [GHz]", color=INK_2)
ax.set_ylabel("HFSS dataset used [GHz]", color=INK_2)
ax.set_title("a. Each photon gets the grid point nearest in log frequency",
             fontsize=11, color=INK, loc="left")
ax.legend(handles=[Line2D([], [], marker="o", linestyle="none", markersize=5,
                          color=INK_MUTED, label="broadband 20 K photons"),
                   Line2D([], [], marker="o", linestyle="none", markersize=8,
                          markerfacecolor="none", markeredgecolor=INK,
                          label="fixed-gun probe runs"),
                   Line2D([], [], linestyle="--", color=INK_MUTED,
                          label="log-space midpoints")],
          frameon=False, fontsize=8, labelcolor=INK_2, loc="upper left")

# --- (b) transmittance, observed vs predicted --------------------------------
ax = axes[0, 1]
grids = hfss.discover_frequencies(CRACK1, wg_dir)
stems = dict(grids)
rows = []
for rid, g in gun.groupby("run_id"):
    f_sel = float(g["hfss_freq_GHz"].iloc[0])
    n = len(g)
    k = int((g["status"] == "BBRDiffractionTransmit").sum())
    ds = hfss.load_dataset(stems[f_sel], wg_dir)
    inc = hfss.fold_incidence(g[["px_pre", "py_pre", "pz_pre"]].to_numpy(float)[0])
    t_pred, _ = hfss.random_polarization_mixture(ds[hfss.nearest_key(
        ds, inc.phi_deg, inc.theta_deg)])
    rows.append((f_sel, k / n, np.sqrt((k / n) * (1 - k / n) / n), t_pred))
obs = np.array(rows)
jitter = {f: 0 for f in GRID}
for f_sel, t_obs, sig, t_pred in rows:
    x = GRID.index(f_sel) + (jitter[f_sel] - 1) * 0.11
    jitter[f_sel] += 1
    ax.errorbar(x, t_obs, yerr=3 * sig, fmt="o", markersize=6,
                color=COLOR[f_sel], ecolor=COLOR[f_sel], elinewidth=1.2,
                capsize=0, zorder=3)
for i, f in enumerate(GRID):
    t_pred = [r[3] for r in rows if r[0] == f][0]
    ax.hlines(t_pred, i - 0.34, i + 0.34, color=INK, linewidth=1.6, zorder=2)
    ax.annotate(f"{t_pred:.3f}", (i + 0.36, t_pred), fontsize=8, color=INK_2,
                va="center")
ax.set_xticks(range(len(GRID)))
ax.set_xticklabels([f"{f:g}" for f in GRID])
ax.set_xlim(-0.6, len(GRID) - 0.25)
ax.set_xlabel("HFSS dataset used [GHz]", color=INK_2)
ax.set_ylabel("crack transmittance", color=INK_2)
ax.set_title("b. Observed transmittance matches the dataset that was selected",
             fontsize=11, color=INK, loc="left")
ax.legend(handles=[Line2D([], [], marker="o", linestyle="none", markersize=6,
                          color=INK_MUTED, label="observed, 3σ bars (one per run)"),
                   Line2D([], [], color=INK, linewidth=1.6,
                          label="predicted from the mock CSVs")],
          frameon=False, fontsize=8, labelcolor=INK_2, loc="upper left")

# --- (c) broadband spectrum, coloured by dataset -----------------------------
ax = axes[1, 0]
bins = np.logspace(np.log10(max(planck["nu_GHz"].min(), 5.0)),
                   np.log10(planck["nu_GHz"].max()), 90)
ax.hist([planck[planck["hfss_freq_GHz"] == f]["nu_GHz"] for f in GRID],
        bins=bins, stacked=True, color=[COLOR[f] for f in GRID],
        label=[f"{f:g} GHz" for f in GRID], linewidth=0)
for m in MIDS:
    ax.axvline(m, color=INK_MUTED, linestyle="--", linewidth=0.9)
ax.set_xscale("log")
ax.set_xlabel("photon frequency [GHz]", color=INK_2)
ax.set_ylabel("photons entering crack1", color=INK_2)
ax.set_title("c. A 20 K Planck run spreads across the whole grid",
             fontsize=11, color=INK, loc="left")
leg = ax.legend(frameon=False, fontsize=8, labelcolor=INK_2, title="dataset used",
                loc="upper left")
leg.get_title().set_color(INK_2); leg.get_title().set_fontsize(8)

# --- (d) far-field signature -------------------------------------------------
ax = axes[1, 1]
# Fixed-gun runs only. Those are at normal incidence, where the fold leaves the
# axes unmirrored and a far-field Theta <= 90 deg maps to k_z >= 0 in the world
# frame. The broadband run arrives at all azimuths, and the quarter-symmetry
# unfold flips phi_hat for half of them, so its k_z carries no signature.
tx = gun[gun["status"] == "BBRDiffractionTransmit"]
fracs = [( (tx[tx["hfss_freq_GHz"] == f]["pz_post"] < 0).mean() if
           (tx["hfss_freq_GHz"] == f).any() else np.nan) for f in GRID]
bars = ax.bar(range(len(GRID)), fracs, width=0.55,
              color=[COLOR[f] for f in GRID], linewidth=0)
for i, (f, v) in enumerate(zip(GRID, fracs)):
    ax.annotate(f"{v:.3f}", (i, v), textcoords="offset points", xytext=(0, 4),
                ha="center", fontsize=9, color=INK_2)
    if f in TRUNCATED:
        ax.annotate("far field\ntruncated", (i, 0.055), ha="center", fontsize=8,
                    color=INK_2, va="bottom")
ax.axhline(0.5, color=INK_MUTED, linestyle="--", linewidth=0.9)
ax.annotate("symmetric far field", (len(GRID) - 0.55, 0.5),
            textcoords="offset points", xytext=(0, 5), fontsize=8, color=INK_2)
ax.set_xticks(range(len(GRID)))
ax.set_xticklabels([f"{f:g}" for f in GRID])
ax.set_ylim(0, 0.62)
ax.set_xlabel("HFSS dataset used [GHz]", color=INK_2)
ax.set_ylabel("transmitted photons with $k_z < 0$", color=INK_2)
ax.set_title("d. Exit directions come from that dataset too (gun runs, normal incidence)",
             fontsize=11, color=INK, loc="left")

fig.suptitle("Frequency-keyed HFSS lookup, validated against a mock five-frequency tree",
             fontsize=13, color=INK, x=0.011, ha="left", y=0.985)
fig.tight_layout(rect=(0, 0, 1, 0.965))
fig.savefig(out_path, dpi=150, facecolor=SURFACE, bbox_inches="tight")
print(f"wrote {out_path}")

print(f"\n{'dataset':>8} {'runs':>5} {'T_obs':>18} {'T_pred':>8} {'k_z<0':>7}")
for i, f in enumerate(GRID):
    rs = [r for r in rows if r[0] == f]
    o = ", ".join(f"{r[1]:.3f}" for r in rs)
    print(f"{f:>8g} {len(rs):>5} {o:>18} {rs[0][3]:>8.3f} {fracs[i]:>7.3f}")
