"""
Section 1: Replicate Joey Bernhardt's (JB) direct P* estimation pipeline
in Python using her raw phosphate RFU data.

JB method (12c_phosphate_growth_direct.R):
  1. ID exponential growth phase per well_plate (sliding log-linear windows 3-7 pts)
  2. N0_mean = mean first-RFU across all well_plates within each population
  3. Pool all 4 replicates per population (group by ancestor_id, treatment, population)
  4. Fit Monod directly to raw RFU timeseries:
       RFU ~ N0_mean * exp(umax * P/(ks+P) * days)
  5. R* = ks * m / (umax - m),  m = 0.56 day^-1
  6. Change-from-ancestor = R* - mean(ancestral R*)
  7. Compare to published change-pstar-monod-boot-0.56.csv  ->  R^2

Output files (all in JB-v-JL/):
  JB_replicated_rstars.csv
  JB_replicated_change_pstar.csv
  section1_comparison_boot.csv
  section1_comparison_pt.csv
  section1_fig1.png
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# NLS fitting (numpy only -- no scipy required)
# ---------------------------------------------------------------------------

def monod_pred(days, P, N0, umax, ks):
    return N0 * np.exp((umax * P / (ks + P)) * days)


def fit_monod(days, P, N0, RFU,
              umax_bounds=(0, 3), ks_bounds=(0, 100),
              n_ks=300, max_iter=300):
    """
    Fit RFU ~ N0*exp(umax*P/(ks+P)*days).

    Phase 1 : profile over ks grid with linear solve for umax in log-space.
    Phase 2 : Gauss-Newton / LM refinement in original RFU space.
    Returns (umax, ks).
    """
    days = np.asarray(days, float)
    P    = np.asarray(P,    float)
    N0   = np.asarray(N0,   float)
    y    = np.asarray(RFU,  float)

    # --- Phase 1: profile over ks ---
    safe = (y > 0) & (N0 > 0)
    log_y = np.log(y[safe]) - np.log(N0[safe])

    ks_lo = max(ks_bounds[0], 1e-6)
    ks_grid = np.linspace(ks_lo, ks_bounds[1], n_ks)
    best_ssr, best_ks, best_umax = np.inf, ks_grid[0], 1.0

    for ks in ks_grid:
        z = P[safe] / (ks + P[safe]) * days[safe]
        denom = np.dot(z, z)
        if denom < 1e-12:
            continue
        umax = np.dot(z, log_y) / denom
        umax = float(np.clip(umax, umax_bounds[0] + 1e-9, umax_bounds[1]))
        ssr = float(np.sum((log_y - umax * z) ** 2))
        if ssr < best_ssr:
            best_ssr, best_ks, best_umax = ssr, float(ks), umax

    # --- Phase 2: Gauss-Newton / LM refinement ---
    p   = np.array([best_umax, best_ks], dtype=float)
    lo  = np.array([umax_bounds[0], ks_bounds[0]], dtype=float)
    hi  = np.array([umax_bounds[1], ks_bounds[1]], dtype=float)
    lam = 1e-3

    for _ in range(max_iter):
        pred  = monod_pred(days, P, N0, p[0], p[1])
        resid = y - pred
        z_frac = P / (p[1] + P)
        dfu = pred * z_frac * days
        dfk = -pred * p[0] * days * P / (p[1] + P) ** 2
        J   = np.column_stack([dfu, dfk])

        JTJ  = J.T @ J
        JTr  = J.T @ resid
        reg  = lam * np.diag(np.diag(JTJ) + 1e-10)
        try:
            dp = np.linalg.solve(JTJ + reg, JTr)
        except np.linalg.LinAlgError:
            break

        p_new = np.clip(p + dp, lo, hi)
        new_ssr = float(np.sum((y - monod_pred(days, P, N0, p_new[0], p_new[1])) ** 2))
        old_ssr = float(np.sum(resid ** 2))

        if new_ssr < old_ssr:
            p   = p_new
            lam = max(lam / 10, 1e-10)
        else:
            lam = min(lam * 10, 1e10)

        if np.linalg.norm(dp) < 1e-9:
            break

    return float(p[0]), float(p[1])


def pearsonr(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    xm = x - x.mean()
    ym = y - y.mean()
    denom = np.sqrt(np.dot(xm, xm) * np.dot(ym, ym))
    return float(np.dot(xm, ym) / denom) if denom > 0 else 0.0


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE    = "/sessions/vibrant-pensive-wozniak/mnt/Chlamy_37_pops_pheno_constraint"
JB_DATA = BASE + "/JB-v-JL/chlamee-r-star/data-processed"
OUT     = BASE + "/JB-v-JL"

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("Loading RFU timeseries ...")
rfu = pd.read_csv(JB_DATA + "/phosphate-rstar-rfus-time.csv")

# Treatment / ancestor mapping from JB's published point estimates
treat_map = pd.read_csv(JB_DATA + "/phosphate-rstars-direct.csv")[
    ['ancestor_id', 'treatment', 'population']
].copy()
treat_map['population'] = treat_map['population'].astype(str)

# JB published results
jb_pub  = pd.read_csv(JB_DATA + "/change-pstar-monod-boot-0.56.csv")   # bootstrap
jb_pt   = pd.read_csv(JB_DATA + "/phosphate-rstars-direct.csv")         # point (m=0.1)

# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------
rfu = rfu[rfu['population'] != 'COMBO'].copy()
rfu = rfu[rfu['well_plate'] != 'D11_30'].copy()     # JB's outlier removal
rfu = rfu.dropna(subset=['RFU'])
rfu['population'] = rfu['population'].astype(str)
rfu = rfu.sort_values(['well_plate', 'days'])

# N0 = first RFU per well_plate
n0_per_well = (rfu.groupby('well_plate')
                   .apply(lambda g: g.iloc[0]['RFU'])
                   .reset_index(name='N0'))
rfu = rfu.merge(n0_per_well, on='well_plate', how='left')

# N0_mean = mean N0 across well_plates within the same population
n0_mean = (rfu.groupby('population')['N0']
               .mean()
               .reset_index(name='N0_mean'))
rfu = rfu.merge(n0_mean, on='population', how='left')

# ---------------------------------------------------------------------------
# Identify exponential phase (sliding log-linear windows 3-7 points)
# ---------------------------------------------------------------------------
print("Identifying exponential windows ...")

def best_window(group):
    group = group.sort_values('days').reset_index(drop=True)
    best_slope, best_n = -np.inf, 3
    for n in range(3, 8):
        sub = group.head(n)
        if sub['RFU'].min() <= 0 or len(sub) < 2:
            continue
        x = sub['days'].values
        y = np.log(sub['RFU'].values)
        slope = np.polyfit(x, y, 1)[0]
        if slope > best_slope:
            best_slope, best_n = slope, n
    return best_n

window_df = (rfu.groupby('well_plate')
                 .apply(best_window)
                 .reset_index(name='n_pts'))

rfu['time_point'] = rfu.groupby('well_plate').cumcount() + 1
rfu = rfu.merge(window_df, on='well_plate', how='left')
rfu_exp = rfu[rfu['time_point'] <= rfu['n_pts']].copy()

# Merge treatment info
rfu_exp = rfu_exp.merge(treat_map, on='population', how='left')
print("  Rows with treatment info: %d / %d" % (rfu_exp['ancestor_id'].notna().sum(), len(rfu_exp)))

# ---------------------------------------------------------------------------
# Direct Monod fitting (replicates pooled per population -- JB approach)
# ---------------------------------------------------------------------------
MORTALITY = 0.56

print("Fitting Monod curves ...")
results = []

for (ancestor_id, treatment, population), grp in (
        rfu_exp.dropna(subset=['ancestor_id'])
               .groupby(['ancestor_id', 'treatment', 'population'])):

    g = grp.dropna(subset=['RFU', 'N0_mean', 'phosphate_concentration', 'days'])
    g = g[(g['RFU'] > 0) & (g['N0_mean'] > 0)]

    if len(g) < 5:
        print("  Skip %s (%s): only %d pts" % (population, treatment, len(g)))
        continue

    try:
        umax, ks = fit_monod(
            g['days'].values,
            g['phosphate_concentration'].values,
            g['N0_mean'].values,
            g['RFU'].values
        )
        rstar = ks * MORTALITY / (umax - MORTALITY) if umax > MORTALITY else np.nan
        results.append({
            'ancestor_id' : ancestor_id,
            'treatment'   : treatment,
            'population'  : population,
            'umax'        : umax,
            'ks'          : ks,
            'rstar'       : rstar,
        })
    except Exception as e:
        print("  Fit failed %s (%s): %s" % (population, treatment, e))

results_df = pd.DataFrame(results)
print("  Successfully fit: %d populations" % len(results_df))

# ---------------------------------------------------------------------------
# Change from ancestral mean R*
# ---------------------------------------------------------------------------
anc_rstars = (results_df[results_df['treatment'] == 'none']
              .groupby('ancestor_id')['rstar']
              .mean()
              .reset_index(name='anc_rstar_mean'))

results_df = results_df.merge(anc_rstars, on='ancestor_id', how='left')
results_df['change_rstar'] = results_df['rstar'] - results_df['anc_rstar_mean']
results_df.to_csv(OUT + "/JB_replicated_rstars.csv", index=False)

# ---------------------------------------------------------------------------
# Comparison against JB published bootstrap means
# ---------------------------------------------------------------------------
print("\nComparing to JB published estimates ...")

jb_pub2 = jb_pub.copy()
jb_pub2['population'] = jb_pub2['population'].astype(str)
jb_pub2['treatment_code'] = jb_pub2['treatment'].replace({'Ancestors': 'none'})

comp_boot = results_df.merge(
    jb_pub2[['population', 'ancestor_id', 'treatment_code',
             'change_rstar_mean', 'change_rstar_lower', 'change_rstar_upper']],
    left_on  =['population', 'ancestor_id', 'treatment'],
    right_on =['population', 'ancestor_id', 'treatment_code'],
    how='inner'
)
comp_boot.to_csv(OUT + "/section1_comparison_boot.csv", index=False)

# Compare against JB point estimates re-computed with m=0.56
jb_pt2 = jb_pt.copy()
jb_pt2['rstar_056'] = jb_pt2['ks'] * 0.56 / (jb_pt2['umax'] - 0.56)
jb_anc = (jb_pt2[jb_pt2['treatment'] == 'none']
          .groupby('ancestor_id')['rstar_056']
          .mean()
          .reset_index(name='anc_rstar_jb'))
jb_pt2 = jb_pt2.merge(jb_anc, on='ancestor_id', how='left')
jb_pt2['change_rstar_jb_pt'] = jb_pt2['rstar_056'] - jb_pt2['anc_rstar_jb']
jb_pt2['population'] = jb_pt2['population'].astype(str)

comp_pt = results_df.merge(
    jb_pt2[['population', 'ancestor_id', 'treatment',
             'umax', 'ks', 'rstar_056', 'change_rstar_jb_pt']],
    on=['population', 'ancestor_id', 'treatment'],
    suffixes=('_rep', '_jb')
)
comp_pt.to_csv(OUT + "/section1_comparison_pt.csv", index=False)

# R^2 values
valid_b = comp_boot.dropna(subset=['change_rstar', 'change_rstar_mean'])
r_b  = pearsonr(valid_b['change_rstar'], valid_b['change_rstar_mean'])
r2_b = r_b ** 2
print("  R^2 (rep vs JB bootstrap mean): %.4f  (r=%.4f, n=%d)" % (r2_b, r_b, len(valid_b)))

valid_p = comp_pt.dropna(subset=['change_rstar', 'change_rstar_jb_pt'])
r_p  = pearsonr(valid_p['change_rstar'], valid_p['change_rstar_jb_pt'])
r2_p = r_p ** 2
print("  R^2 (rep vs JB point est m=0.56): %.4f  (n=%d)" % (r2_p, len(valid_p)))

# umax check
valid_u  = comp_pt.dropna(subset=['umax_rep', 'umax_jb'])
r_u  = pearsonr(valid_u['umax_jb'], valid_u['umax_rep'])
r2_u = r_u ** 2

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
TCOLS = {
    'B': '#f9a729', 'BS': '#97cfd0', 'C': '#00a2b3',
    'L': '#f1788d', 'N':  '#cf3e53', 'P': '#b9ca5d',
    'S': '#6b6b6b', 'none': 'black',  'Ancestors': 'black'
}

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# Panel A -- change in R* (replicated vs published)
ax = axes[0]
for trt, grp in valid_b.groupby('treatment'):
    ax.scatter(grp['change_rstar_mean'], grp['change_rstar'],
               color=TCOLS.get(trt, 'grey'), label=trt,
               s=60, alpha=0.85, edgecolors='k', linewidths=0.5)

lim = max(abs(valid_b['change_rstar']).max(),
          abs(valid_b['change_rstar_mean']).max()) * 1.15
ax.plot([-lim, lim], [-lim, lim], 'k--', lw=1, alpha=0.5)
ax.axhline(0, color='grey', lw=0.5, ls=':')
ax.axvline(0, color='grey', lw=0.5, ls=':')
ax.set_xlabel("JB published delta-P* (bootstrap mean, uM)", fontsize=11)
ax.set_ylabel("JB replicated delta-P* (Python NLS, uM)", fontsize=11)
ax.set_title("Section 1: Replication check\nR^2 = %.3f, n = %d" % (r2_b, len(valid_b)), fontsize=11)
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
handles = [plt.scatter([], [], color=TCOLS.get(t, 'grey'), s=40, label=t,
                       edgecolors='k', linewidths=0.4)
           for t in sorted(valid_b['treatment'].unique()) if t != 'none']
ax.legend(handles=handles, title='Treatment', fontsize=8, title_fontsize=8,
          framealpha=0.7, loc='upper left')

# Panel B -- umax replication
ax2 = axes[1]
for trt, grp in valid_u.groupby('treatment'):
    ax2.scatter(grp['umax_jb'], grp['umax_rep'],
                color=TCOLS.get(trt, 'grey'),
                s=60, alpha=0.85, edgecolors='k', linewidths=0.5)

lim2 = max(valid_u['umax_jb'].max(), valid_u['umax_rep'].max()) * 1.1
ax2.plot([0, lim2], [0, lim2], 'k--', lw=1, alpha=0.5)
ax2.set_xlabel("JB published umax (per day)", fontsize=11)
ax2.set_ylabel("JB replicated umax (per day)", fontsize=11)
ax2.set_title("umax replication\nR^2 = %.3f, n = %d" % (r2_u, len(valid_u)), fontsize=11)

plt.tight_layout()
fig.savefig(OUT + "/section1_fig1.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: " + OUT + "/section1_fig1.png")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print("\n" +"="*60)
print("SECTION 1 SUMMARY")
print("="*60)
print("  Populations fit    : %d" % len(results_df))
print("  R^2 vs JB boot mean: %.4f" % r2_b)
print("  R^2 vs JB pt m=0.56: %.4f" % r2_p)
if r2_b > 0.95:
    print("  --> Excellent replication")
elif r2_b > 0.85:
    print("  --> Good replication")
else:
    print("  --> Moderate -- check individual fits")
print()
print("First 10 replicated R* values:")
print(results_df[['ancestor_id','treatment','population',
                   'umax','ks','rstar','change_rstar']].head(10).to_string(index=False))
))
