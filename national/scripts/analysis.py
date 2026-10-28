"""Analyses and figures for the iNaturalist mosquito / WNV project.

Reads data/processed/*.csv, writes figures/*.png and results/results.json + results/tables.md
"""
import json, os
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P = os.path.join(ROOT, "data", "processed")
FIG = os.path.join(ROOT, "figures"); os.makedirs(FIG, exist_ok=True)
RES = os.path.join(ROOT, "results"); os.makedirs(RES, exist_ok=True)
YEARS = list(range(2016, 2026))
R = {}

# plain, consistent styling
BLUE, ORANGE, GREY = "#2b6cb0", "#c05621", "#8a8a8a"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})

obs = pd.read_csv(os.path.join(P, "mosquito_obs.csv"), parse_dates=["date"])
nat = pd.read_csv(os.path.join(P, "national_month.csv"))
sy = pd.read_csv(os.path.join(P, "state_year.csv"), dtype={"state_fips": str})
county = pd.read_csv(os.path.join(P, "county.csv"), dtype={"fips": str})
obs = obs[obs["year"].between(2016, 2025)]
nat = nat[nat["year"].between(2016, 2025)].sort_values(["year", "month"]).reset_index(drop=True)


def share(num, den):
    """log of a smoothed share, so zero counts don't blow up"""
    return np.log((num + 0.5) / (den + 1.0))


# =================================================================== 1. descriptives
d = {}
d["n_obs"] = int(len(obs))
d["n_observers"] = int(obs["user"].nunique())
d["pct_culex"] = round(100 * obs["culex"].mean(), 1)
d["pct_family_only"] = round(100 * (obs["rank"] == "family").mean(), 1)
d["pct_research_grade"] = round(100 * (obs["qg"] == "research").mean(), 1)
d["pct_with_county"] = round(100 * obs["fips"].notna().mean(), 1)
per_user = obs.groupby("user").size().sort_values(ascending=False)
d["pct_obs_top1pct_users"] = round(100 * per_user.iloc[:max(1, len(per_user) // 100)].sum() / len(obs), 1)
d["pct_obs_top10_users"] = round(100 * per_user.iloc[:10].sum() / len(obs), 1)
d["pct_users_single_obs"] = round(100 * (per_user == 1).mean(), 1)
yr = obs.groupby("year").agg(mosq_obs=("id", "size"), observers=("user", "nunique"))
yr = yr.join(nat.groupby("year")[["insect_obs", "wnv_cases", "wnv_neuro"]].sum())
d["by_year"] = yr.reset_index().to_dict(orient="records")
R["descriptives"] = d

fig, ax = plt.subplots(figsize=(6.5, 3.2))
ax.plot(yr.index, yr["mosq_obs"], "-o", color=BLUE, label="iNaturalist mosquito observations")
ax.set_ylabel("mosquito observations", color=BLUE)
ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
ax2.plot(yr.index, yr["wnv_neuro"], "-s", color=ORANGE, label="WNV neuroinvasive cases (CDC)")
ax2.set_ylabel("neuroinvasive WNV cases", color=ORANGE)
ax.set_title("Crowd activity grows every year; WNV does not")
ax.set_xticks(YEARS)
fig.savefig(os.path.join(FIG, "fig1_growth.png")); plt.close(fig)

# =================================================================== 2. national: seasonal timing and year-to-year
nat["mosq_share"] = (nat["mosq_obs"] + 0.5) / (nat["insect_obs"] + 1)
# seasonal profile = share of all 2016-2025 records that fall in each calendar month (pooled counts;
# dominated by recent years, when volumes are large enough to give a stable shape)
cols = ["mosq_obs", "insect_obs", "cranefly_obs", "wnv_neuro"]
prof = nat.groupby("month")[cols].sum()
prof = prof / prof.sum()
peak = {c: int(prof[c].idxmax()) for c in cols}
cent = {c: round(float((prof[c] * prof.index).sum() / prof[c].sum()), 2) for c in cols}
lagfit = {}
for k in range(0, 4):     # mosquito profile shifted k months later vs. WNV profile
    r, _ = stats.pearsonr(np.roll(prof["mosq_obs"].values, k), prof["wnv_neuro"].values)
    lagfit[k] = round(r, 3)
R["seasonal"] = {"peak_month": peak, "center_of_mass_month": cent, "profile_corr_by_shift": lagfit}
# Culex-only timing, from the observation table
cx_prof = obs[obs["culex"]].groupby("month").size().reindex(range(1, 13), fill_value=0)
cx_prof = cx_prof / cx_prof.sum()
R["seasonal"]["monthly_pct"] = {"mosquito": (100 * prof["mosq_obs"]).round(1).tolist(),
                                "culex": (100 * cx_prof).round(1).tolist(),
                                "insect": (100 * prof["insect_obs"]).round(1).tolist(),
                                "wnv_neuro": (100 * prof["wnv_neuro"]).round(1).tolist()}
R["seasonal"]["culex_center_of_mass_month"] = round(float((cx_prof * cx_prof.index).sum()), 2)
R["seasonal"]["culex_peak_month"] = int(cx_prof.idxmax())
gpath = os.path.join(ROOT, "data", "raw", "inat", "genus_counts.json")
if os.path.exists(gpath):
    g = json.load(open(gpath))
    R["genera_pct_of_mosquito_obs"] = {k: round(100 * v["obs_2016_2025"] / len(obs), 1) for k, v in g.items()}

# month-level anomalies, May-Oct. x = log mosquito share of insect obs, y = log(1+neuro cases),
# each minus its calendar-month mean over 2016-2025.
ly = np.log1p(nat["wnv_neuro"])
lx = np.log(nat["mosq_share"])
nat["x_anom"] = lx - lx.groupby(nat["month"]).transform("mean")
nat["y_anom"] = ly - ly.groupby(nat["month"]).transform("mean")
rng = np.random.default_rng(0)


def paired(df, k, within_year):
    xs, ys = [], []
    for y0 in YEARS:
        a = df[df["year"] == y0].set_index("month")
        xa, ya = a["x_anom"], a["y_anom"]
        if within_year:      # remove the year's average level -> only timing within the season is left
            xa = xa - xa.loc[5:10].mean(); ya = ya - ya.loc[5:10].mean()
        for m_ in range(5, 11):
            if m_ + k <= 12:
                xs.append(xa.loc[m_]); ys.append(ya.loc[m_ + k])
    return np.array(xs), np.array(ys)


def perm_p(k, within_year, r_obs, n=5000):
    """shuffle which year's mosquito series is paired with which year's WNV series"""
    hits = 0
    x_by_year = {y0: nat[nat["year"] == y0].set_index("month")["x_anom"] for y0 in YEARS}
    for _ in range(n):
        perm = dict(zip(YEARS, rng.permutation(YEARS)))
        df = nat.copy()
        df["x_anom"] = [x_by_year[perm[y0]].loc[m_] for y0, m_ in zip(df["year"], df["month"])]
        xs, ys = paired(df, k, within_year)
        hits += abs(stats.pearsonr(xs, ys)[0]) >= abs(r_obs)
    return (hits + 1) / (n + 1)


an = {}
for within in (False, True):
    for k in range(0, 3):
        xs, ys = paired(nat, k, within)
        r, p = stats.pearsonr(xs, ys)
        an[f"{'within_year' if within else 'total'}_lag{k}"] = {
            "r": round(r, 3), "p_naive": round(p, 4), "p_year_permutation": round(perm_p(k, within, r, 2000), 4),
            "n": len(xs)}
R["monthly_anomaly"] = an

# year level: early-season (Apr-Jul) national mosquito share vs. that year's neuroinvasive cases
e = nat[nat["month"].between(4, 7)].groupby("year")[["mosq_obs", "insect_obs", "cranefly_obs"]].sum()
e = e.join(nat.groupby("year")["wnv_neuro"].sum())
e["share"] = 1000 * e["mosq_obs"] / e["insect_obs"]; e["crane"] = 1000 * e["cranefly_obs"] / e["insect_obs"]
yl = {}
for c in ["share", "crane"]:
    r, p = stats.spearmanr(e[c], e["wnv_neuro"]); yl[c] = {"spearman_rho": round(r, 3), "p": round(p, 4), "n": len(e)}
# concurrent window for comparison: Jul-Sep mosquito share (same months as most onsets)
c3 = nat[nat["month"].between(7, 9)].groupby("year")[["mosq_obs", "insect_obs"]].sum()
c3 = (1000 * c3["mosq_obs"] / c3["insect_obs"]).rename("share_JulSep")
e = e.join(c3)
r, p = stats.spearmanr(e["share_JulSep"], e["wnv_neuro"]); yl["share_JulSep"] = {"spearman_rho": round(r, 3), "p": round(p, 4), "n": len(e)}
R["year_level"] = yl
R["year_level_table"] = e.round(3).reset_index().to_dict(orient="records")

fig, axes = plt.subplots(1, 2, figsize=(9, 3.2), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axes[0]
ax.plot(prof.index, prof["insect_obs"], "-", color=GREY, lw=1, label="all insect observations")
ax.plot(prof.index, prof["mosq_obs"], "-o", ms=4, color=BLUE, label="mosquito observations")
ax.plot(cx_prof.index, cx_prof.values, ":", lw=1.2, color=BLUE, label="Culex observations only")
ax.plot(prof.index, prof["wnv_neuro"], "-s", ms=4, color=ORANGE, label="WNV neuroinvasive onsets")
ax.set_xticks(range(1, 13)); ax.set_xticklabels(list("JFMAMJJASOND"))
ax.set_ylabel("share of 2016-2025 total"); ax.set_title("(a) Seasonal timing, 2016-2025 pooled", fontsize=10)
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax = axes[1]
z = lambda v: (v - v.mean()) / v.std()
ax.plot(e.index, z(e["share"]), "-o", ms=4, color=BLUE, label="Apr-Jul mosquito share of insect obs")
ax.plot(e.index, z(e["wnv_neuro"]), "-s", ms=4, color=ORANGE, label="neuroinvasive WNV cases")
ax.axhline(0, color=GREY, lw=.5)
ax.set_xticks(YEARS[::3]); ax.set_ylabel("z-score")
ax.set_title(f"(b) Year to year (Spearman ρ={yl['share']['spearman_rho']:.2f}, n=10)", fontsize=10)
ax.legend(frameon=False, fontsize=7.5, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig2_national.png")); plt.close(fig)

# =================================================================== 3. state-year panel (main test)
for w in ["early", "full"]:
    sy[f"x_share_{w}"] = share(sy[f"mosq_obs_{w}"], sy[f"insect_obs_{w}"])
    sy[f"x_obsr_{w}"] = share(sy[f"mosq_observers_{w}"], sy[f"insect_obs_{w}"])
    sy[f"x_culex_{w}"] = share(sy[f"culex_obs_{w}"], sy[f"insect_obs_{w}"])
    sy[f"x_crane_{w}"] = share(sy[f"cranefly_obs_{w}"], sy[f"insect_obs_{w}"])
    sy[f"x_count_{w}"] = np.log1p(sy[f"mosq_obs_{w}"])
    sy[f"x_obscount_{w}"] = np.log1p(sy[f"mosq_observers_{w}"])
    sy[f"x_insect_{w}"] = np.log1p(sy[f"insect_obs_{w}"])
sy["y"] = np.log1p(sy["neuro"])
sy["y_all"] = np.log1p(sy["cases"])

MIN_INSECT = 200          # early-season insect observations needed for a stable share
panel = sy[sy["insect_obs_early"] >= MIN_INSECT].copy()
R["panel_n"] = {"state_years": int(len(panel)), "states": int(panel["state"].nunique()),
                "dropped_low_activity": int(len(sy) - len(panel))}


def fe(formula, data):
    m = smf.ols(formula + " + C(state) + C(year)", data=data).fit(
        cov_type="cluster", cov_kwds={"groups": data["state"]})
    return m

signals = {
    "x_share_early": "Mosquito share of insect obs (Apr-Jul)",
    "x_obsr_early": "Mosquito observers per insect obs (Apr-Jul)",
    "x_culex_early": "Culex share of insect obs (Apr-Jul)",
    "x_count_early": "Raw mosquito obs count (Apr-Jul)",
    "x_obscount_early": "Unique mosquito observers (Apr-Jul)",
    "x_crane_early": "PLACEBO: crane-fly share of insect obs (Apr-Jul)",
    "x_share_full": "Mosquito share of insect obs (full year, concurrent)",
}
tab = []
for x, label in signals.items():
    m = fe(f"y ~ {x}", panel)
    # within-R2: variance explained beyond state + year effects
    base = smf.ols("y ~ C(state) + C(year)", data=panel).fit()
    within_r2 = 1 - m.ssr / base.ssr
    tab.append({"signal": label, "var": x, "beta": round(m.params[x], 3), "se": round(m.bse[x], 3),
                "p": round(m.pvalues[x], 4), "within_R2": round(within_r2, 4), "n": int(m.nobs)})
R["panel_fe"] = tab

# robustness: control for overall insect activity, and all-disease cases as outcome
m1 = fe("y ~ x_count_early + x_insect_early", panel)
m2 = fe("y_all ~ x_share_early", panel)
R["panel_robust"] = {
    "count_with_insect_control": {"beta_count": round(m1.params["x_count_early"], 3),
                                  "p": round(m1.pvalues["x_count_early"], 4),
                                  "beta_insect": round(m1.params["x_insect_early"], 3)},
    "all_cases_outcome_share": {"beta": round(m2.params["x_share_early"], 3), "p": round(m2.pvalues["x_share_early"], 4)},
}
# iNaturalist volume was small before 2020; re-run the main signals on 2020-2025 only
late = panel[panel["year"] >= 2020]
R["panel_2020_2025"] = {}
for x in ["x_share_early", "x_obsr_early", "x_culex_early", "x_crane_early", "x_share_full"]:
    m = fe(f"y ~ {x}", late)
    R["panel_2020_2025"][x] = {"beta": round(m.params[x], 3), "se": round(m.bse[x], 3),
                               "p": round(m.pvalues[x], 4), "n": int(m.nobs)}
# May-Jun only: a cleaner lead, since some July onsets overlap the Apr-Jul window
mj = sm_mj = None
try:
    smo = pd.read_csv(os.path.join(P, "state_month.csv"), dtype={"state_fips": str})
    smo = smo[smo["month"].isin([5, 6])].groupby(["state", "year"])[["mosq_obs", "insect_obs"]].sum().reset_index()
    mj = panel.merge(smo.rename(columns={"mosq_obs": "m_mj", "insect_obs": "i_mj"}), on=["state", "year"])
    mj["x_share_mj"] = share(mj["m_mj"], mj["i_mj"])
    m = fe("y ~ x_share_mj", mj)
    R["panel_MayJun"] = {"beta": round(m.params["x_share_mj"], 3), "se": round(m.bse["x_share_mj"], 3),
                         "p": round(m.pvalues["x_share_mj"], 4), "n": int(m.nobs)}
except Exception as e:
    R["panel_MayJun"] = {"error": str(e)}

# naive pooled correlations (the trap: both scale with state size and iNat growth)
naive = {}
for x in ["x_count_early", "x_share_early"]:
    r, p = stats.spearmanr(panel[x], panel["y"])
    naive[x] = {"spearman_rho": round(r, 3), "p": round(p, 4)}
R["panel_naive"] = naive

# negative-binomial version on raw counts (state + year FE), share signal
import statsmodels.api as sm
try:
    nb = smf.negativebinomial("neuro ~ x_share_early + C(state) + C(year)", data=panel).fit(disp=0, maxiter=500)
    alpha = nb.params["alpha"]
    glm = smf.glm("neuro ~ x_share_early + C(state) + C(year)", data=panel,
                  family=sm.families.NegativeBinomial(alpha=alpha)).fit(
        cov_type="cluster", cov_kwds={"groups": pd.factorize(panel["state"])[0]})
    R["panel_nb"] = {"beta": round(nb.params["x_share_early"], 3), "p_model_se": round(nb.pvalues["x_share_early"], 4),
                     "beta_glm": round(glm.params["x_share_early"], 3), "p_cluster_se": round(glm.pvalues["x_share_early"], 4),
                     "alpha": round(alpha, 3)}
except Exception as e:
    R["panel_nb"] = {"error": str(e)}

# without year effects: lets nationally bad years count as signal (state effects + linear trend only)
m3 = smf.ols("y ~ x_share_early + C(state) + year", data=panel).fit(cov_type="cluster", cov_kwds={"groups": panel["state"]})
R["panel_no_year_fe"] = {"beta": round(m3.params["x_share_early"], 3), "se": round(m3.bse["x_share_early"], 3),
                         "p": round(m3.pvalues["x_share_early"], 4)}

# per-state within-state rank correlations (sign test)
per = []
for st, g in panel.groupby("state"):
    if len(g) >= 7 and g["neuro"].sum() >= 10:
        r, _ = stats.spearmanr(g["x_share_early"], g["y"])
        per.append((st, len(g), r))
per = pd.DataFrame(per, columns=["state", "years", "rho"]).dropna()
pos = int((per["rho"] > 0).sum())
R["per_state"] = {"n_states": int(len(per)), "n_positive": pos, "median_rho": round(per["rho"].median(), 3),
                  "sign_test_p": round(stats.binomtest(pos, len(per), 0.5).pvalue, 4),
                  "table": per.round(3).to_dict(orient="records")}

# out-of-sample: leave-one-year-out forecast of each state's season
def loyo(xvar):
    errs_b, errs_m, hits = [], [], []
    for y0 in YEARS:
        tr, te = panel[panel["year"] != y0], panel[panel["year"] == y0]
        if te.empty:
            continue
        mu_y = tr.groupby("state")["y"].mean(); mu_x = tr.groupby("state")[xvar].mean()
        trd = tr.assign(yd=tr["y"] - tr["state"].map(mu_y), xd=tr[xvar] - tr["state"].map(mu_x))
        b = np.polyfit(trd["xd"], trd["yd"], 1)[0]
        te = te[te["state"].isin(mu_y.index)]
        base = te["state"].map(mu_y)
        pred = base + b * (te[xvar] - te["state"].map(mu_x))
        errs_b += list((te["y"] - base).abs()); errs_m += list((te["y"] - pred).abs())
    return np.mean(errs_b), np.mean(errs_m), len(errs_b)

oos = {}
for x in ["x_share_early", "x_obsr_early", "x_culex_early", "x_crane_early"]:
    b, mm, n = loyo(x)
    oos[x] = {"mae_baseline": round(b, 4), "mae_with_signal": round(mm, 4),
              "improvement_pct": round(100 * (b - mm) / b, 2), "n": n}
R["out_of_sample"] = oos

# figure: two-way demeaned scatter, mosquito vs placebo
def demean(df, col):
    v = df[col] - df.groupby("state")[col].transform("mean") - df.groupby("year")[col].transform("mean") + df[col].mean()
    return v

fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.3), sharey=True)
yd = demean(panel, "y")
for ax, x, title, col in [(axes[0], "x_share_early", "Mosquitoes (Apr-Jul share)", BLUE),
                          (axes[1], "x_crane_early", "Placebo: crane flies", GREY)]:
    xd = demean(panel, x)
    ax.scatter(xd, yd, s=9, alpha=.45, color=col, edgecolor="none")
    b, a = np.polyfit(xd, yd, 1)
    xs = np.linspace(xd.quantile(.01), xd.quantile(.99), 50)
    ax.plot(xs, a + b * xs, color="black", lw=1)
    row = next(t for t in tab if t["var"] == x)
    ax.set_title(f"{title}\nβ={row['beta']:.2f}, p={row['p']:.3f}", fontsize=9)
    ax.set_xlabel("signal, net of state & year effects")
axes[0].set_ylabel("log(1+neuroinvasive cases),\nnet of state & year effects")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig3_panel.png")); plt.close(fig)

# =================================================================== 4. county cross-section (pooled 2016-2025)
MIN_C = 300
# Connecticut switched to planning-region FIPS codes in 2022, so its old counties are incomplete here
cc = county[(county["insect_obs"] >= MIN_C) & county["pop"].notna() & (county["state"] != "CT")].copy()
cc["x_share"] = share(cc["mosq_obs"], cc["insect_obs"])
cc["x_culex"] = share(cc["culex_obs"], cc["insect_obs"])
cc["inc"] = cc["neuro"] / cc["pop"] * 1e5 / len(YEARS)        # avg annual per 100k
cc["any_neuro"] = (cc["neuro"] > 0).astype(int)
cc["log_pop"] = np.log(cc["pop"]); cc["log_ins"] = np.log(cc["insect_obs"])
cr = {"n_counties": int(len(cc))}
r, p = stats.spearmanr(cc["x_share"], cc["inc"]); cr["spearman_share_incidence"] = {"rho": round(r, 3), "p": round(p, 5)}
r, p = stats.spearmanr(cc["x_culex"], cc["inc"]); cr["spearman_culex_incidence"] = {"rho": round(r, 3), "p": round(p, 5)}
r, p = stats.spearmanr(cc["x_share"], cc["years_any_activity"]); cr["spearman_share_activity_years"] = {"rho": round(r, 3), "p": round(p, 5)}
mc = smf.ols("np.log1p(inc) ~ x_share + log_pop + log_ins + C(state)", data=cc).fit(cov_type="HC1")
cr["ols_within_state"] = {"beta": round(mc.params["x_share"], 3), "p": round(mc.pvalues["x_share"], 4), "n": int(mc.nobs)}
ml = smf.logit("any_neuro ~ x_share + log_pop + log_ins + C(state)", data=cc).fit(disp=0, maxiter=200)
cr["logit_any_case_within_state"] = {"beta": round(ml.params["x_share"], 3), "p": round(ml.pvalues["x_share"], 4)}
R["county"] = cr

cc["bin"] = pd.qcut(cc["x_share"], 10, labels=False, duplicates="drop")
b = cc.groupby("bin").agg(x=("x_share", "median"), inc=("inc", "mean"), n=("fips", "size"),
                          se=("inc", lambda s: s.std() / np.sqrt(len(s))))
fig, ax = plt.subplots(figsize=(5.5, 3.2))
ax.errorbar(np.exp(b["x"]) * 1000, b["inc"], yerr=1.96 * b["se"], fmt="o-", color=BLUE, capsize=2)
ax.set_xscale("log")
ax.set_xticks([0.5, 1, 2, 4, 8]); ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.minorticks_off()
ax.set_xlabel("mosquito observations per 1,000 insect observations (decile medians)")
ax.set_ylabel("avg annual neuroinvasive\nWNV per 100k")
ax.set_title(f"Counties, pooled 2016-2025 (n={len(cc)})", fontsize=10)
fig.savefig(os.path.join(FIG, "fig4_county.png")); plt.close(fig)

# =================================================================== write
json.dump(R, open(os.path.join(RES, "results.json"), "w"), indent=1, default=float)
with open(os.path.join(RES, "tables.md"), "w") as f:
    f.write("| Signal | beta | SE | p | within-R2 | n |\n|---|---|---|---|---|---|\n")
    for t in tab:
        f.write(f"| {t['signal']} | {t['beta']} | {t['se']} | {t['p']} | {t['within_R2']} | {t['n']} |\n")
print(json.dumps(R, indent=1, default=float)[:6000])
