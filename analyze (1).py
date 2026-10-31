"""Are gravitational-wave catalogues a fair map of the universe? (GWTC, O1-O3b)"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

df = pd.read_csv("data/gwtc_confident.csv")
df["mc_det"] = df.mc_src * (1 + df.z)          # what the detector actually "hears"
x, y = np.log10(df.mc_det), np.log10(df.dist_mpc)

# --- Fit log D = a + b log Mc_det, bootstrap the slope -----------------------
b, a, r, p, _ = stats.linregress(x, y)
rho, p_rho = stats.spearmanr(df.mc_det, df.dist_mpc)
rng = np.random.default_rng(42)
boots = []
for _ in range(5000):
    i = rng.integers(0, len(df), len(df))
    boots.append(stats.linregress(x.iloc[i], y.iloc[i])[0])
lo, hi = np.percentile(boots, [2.5, 97.5])

# --- Run-by-run comparison ----------------------------------------------------
order = ["O1/O2", "O3a", "O3b"]
summ = df.groupby("run").agg(n=("name", "size"),
        med_dist=("dist_mpc", "median"), med_mc=("mc_det", "median"),
        med_snr=("snr", "median")).loc[order]
print(summ.round(1)); print(f"\nN={len(df)}  slope={b:.2f} [{lo:.2f},{hi:.2f}]  "
      f"Spearman rho={rho:.2f} (p={p_rho:.1e})  R^2={r**2:.2f}")
# Slope for events with detector-frame Mc < 30 (inspiral-dominated regime)
lowm = df.mc_det < 30
b_lo = stats.linregress(x[lowm], y[lowm])[0]
print(f"slope for Mc_det<30 (n={lowm.sum()}): {b_lo:.2f}")

# --- Figure ------------------------------------------------------------------
cols = {"O1/O2": "#1b6ca8", "O3a": "#e08a1e", "O3b": "#3a9a5b"}
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1.5, 1]})
for run in order:
    d = df[df.run == run]
    ax[0].scatter(d.mc_det, d.dist_mpc, s=d.snr * 6, c=cols[run], alpha=.75,
                  edgecolor="white", lw=.5, label=f"{run} (n={len(d)})")
xx = np.logspace(0, 2.2, 50)
ax[0].plot(xx, 10**a * xx**b, "k-", lw=1.5, label=f"fit: D ∝ M$^{{{b:.2f}}}$")
ref = 10**(y.median()) * (xx / df.mc_det.median())**(5/6)
ax[0].plot(xx, ref, "k--", lw=1, alpha=.6, label="inspiral theory: D ∝ M$^{5/6}$")
ax[0].set(xscale="log", yscale="log", xlabel="Detector-frame chirp mass (M$_\\odot$)",
          ylabel="Luminosity distance (Mpc)",
          title="Heavier systems are heard from farther away")
ax[0].legend(fontsize=8, loc="lower right"); ax[0].grid(alpha=.25, which="both")
ax[0].annotate("GW170817\n(neutron stars)", (df.mc_det[8], 40), (1.8, 250),
               fontsize=8, arrowprops=dict(arrowstyle="-", lw=.6))
ax[0].text(0.02, 0.97, "marker size ∝ SNR", transform=ax[0].transAxes,
           va="top", fontsize=8, color="gray")

for i, run in enumerate(order):
    d = df[df.run == run].dist_mpc
    ax[1].boxplot(d, positions=[i], widths=.5, showfliers=False,
                  medianprops=dict(color="k"))
    ax[1].scatter(np.random.default_rng(i).normal(i, .06, len(d)), d, s=14,
                  c=cols[run], alpha=.7)
ax[1].set(xticks=range(3), xticklabels=order, ylabel="Luminosity distance (Mpc)",
          title="Distance by observing run")
ax[1].grid(alpha=.25, axis="y")
fig.suptitle("GWTC confident events, O1–O3b: the catalogue is a map of sensitivity",
             fontsize=12, y=1.0)
fig.tight_layout(); fig.savefig("figure.png", dpi=160, bbox_inches="tight")
