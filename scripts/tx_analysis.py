"""Texas analysis: do iNaturalist mosquito sightings lead WNV-positive mosquito pools (DSHS)?

Reads   data_tx/tx_panel.csv, data_tx/inat_mosquito_observations_tx.csv
Writes  data_tx/tx_results.json (every number quoted in docs/texas_findings.md),
        figures/tx_seasonal.png, figures/tx_lags.png

Design follows national/scripts/analysis.py:
  - signals are shares of overall iNaturalist activity (all observations; no insect baseline was pulled
    for Texas), on a log scale, over a trailing 4-week window because single county-weeks are mostly empty
  - county + year-week fixed effects, so only differences between counties in the same week count
  - placebo: pair each county's signal with a different year's pools
  - forecast: leave-one-year-out against each county's own history
Run from the repository root.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data_tx")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
rng = np.random.default_rng(0)
R = {}
BLUE, ORANGE, GREY = "#2b6cb0", "#c05621", "#8a8a8a"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})

SEASON = (20, 42)         # MMWR weeks with enough positive pools to measure
WINDOW = 4                # trailing weeks pooled into one signal value
MIN_ALL = 20000           # a county needs this many iNaturalist observations over 2019-2025
LAGS = range(0, 11)
N_PERM = 500

t = pd.read_csv(os.path.join(D, "tx_panel.csv"))
t = t[t["week"] <= 52].sort_values(["county", "year", "week"]).reset_index(drop=True)


def share(num, den):
    return np.log((num + 0.5) / (den + 1.0))


# ----------------------------------------------------------------- signals
g = t.groupby(["county", "year"])
for c in ["all_obs", "mosq_obs", "culex_obs", "mosq_observers"]:
    t[c + "_w"] = g[c].transform(lambda s: s.rolling(WINDOW, min_periods=1).sum())
t["x_share"] = share(t["mosq_obs_w"], t["all_obs_w"])
t["x_obsr"] = share(t["mosq_observers_w"], t["all_obs_w"])
t["x_culex"] = share(t["culex_obs_w"], t["all_obs_w"])
t["y_pools"] = np.log1p(t["new_pools"])
t["y_human"] = np.log1p(t["new_human_total"])

# ----------------------------------------------------------------- sample
tot = t.groupby("county").agg(all_obs=("all_obs", "sum"), mosq=("mosq_obs", "sum"), pools=("new_pools", "sum"))
big = set(tot.index[tot["all_obs"] >= MIN_ALL])
tested = set(tot.index[tot["pools"] > 0])
elig = sorted(big & tested)
R["sample"] = {"counties_total": int(len(tot)), "counties_with_positive_pool": int(len(tested)),
               "counties_ge_min_obs": int(len(big)), "counties_analyzed": int(len(elig)),
               "share_of_pools_in_analyzed": round(float(tot.loc[elig, "pools"].sum() / tot["pools"].sum()), 3),
               "share_of_mosq_obs_in_analyzed": round(float(tot.loc[elig, "mosq"].sum() / tot["mosq"].sum()), 3),
               "min_all_obs": MIN_ALL, "window_weeks": WINDOW, "season_weeks": list(SEASON)}
panel_all = t[t["county"].isin(elig)].copy()

# ----------------------------------------------------------------- descriptives
obs = pd.read_csv(os.path.join(D, "inat_mosquito_observations_tx.csv"), parse_dates=["observed_on"])
per_user = obs.groupby("user_id").size().sort_values(ascending=False)
d = {"mosquito_obs": int(len(obs)), "observers": int(per_user.size),
     "pct_culex": round(100 * obs["is_culex"].mean(), 1),
     "pct_research_grade": round(100 * (obs["quality_grade"] == "research").mean(), 1),
     "pct_top10_users": round(100 * per_user.iloc[:10].sum() / len(obs), 1),
     "pct_users_single_obs": round(100 * (per_user == 1).mean(), 1),
     "pct_county_weeks_with_mosq": round(100 * (t["mosq_obs"] > 0).mean(), 1),
     "pct_county_weeks_with_pool": round(100 * (t["new_pools"].fillna(0) > 0).mean(), 1)}
heavy = per_user[per_user >= 50].index
d["heavy_posters_50plus"] = int(len(heavy))
d["pct_obs_by_heavy_posters"] = round(100 * obs["user_id"].isin(heavy).mean(), 1)
# county-weeks where one person supplied almost everything
cw = obs.assign(week_start=obs["observed_on"].dt.to_period("W-SAT").dt.start_time)
agg = cw.groupby(["county", "week_start"]).agg(n=("id", "size"), top=("user_id", lambda s: s.value_counts().iloc[0]))
agg = agg[agg["n"] >= 5]
d["county_weeks_ge5_sightings"] = int(len(agg))
d["pct_of_those_with_one_user_ge80pct"] = round(100 * float((agg["top"] / agg["n"] >= 0.8).mean()), 1) if len(agg) else None
R["descriptives"] = d

# seasonal timing, pooled 2019-2025 (share of each series' total in each week)
sea = t.groupby("week")[["mosq_obs", "culex_obs", "all_obs", "new_pools", "new_human_total"]].sum()
prof = sea / sea.sum()
cm = lambda s: round(float((s * s.index).sum() / s.sum()), 1)
R["seasonal"] = {"peak_week": {c: int(prof[c].idxmax()) for c in prof}, "center_of_mass_week": {c: cm(sea[c]) for c in sea}}
R["seasonal"]["profile_corr_with_pools"] = {c: round(float(np.corrcoef(prof[c], prof["new_pools"])[0, 1]), 3)
                                            for c in ["mosq_obs", "culex_obs", "all_obs"]}
lagc = {}
for k in range(0, 9):                     # mosquito profile shifted k weeks later
    lagc[k] = round(float(np.corrcoef(np.roll(prof["mosq_obs"].values, k), prof["new_pools"].values)[0, 1]), 3)
R["seasonal"]["profile_corr_by_shift_weeks"] = lagc

fig, ax = plt.subplots(figsize=(6.5, 3.2))
ax.plot(prof.index, prof["all_obs"], color=GREY, lw=1, label="all iNaturalist observations")
ax.plot(prof.index, prof["mosq_obs"], "-o", ms=3, color=BLUE, label="mosquito sightings")
ax.plot(prof.index, prof["culex_obs"], ":", color=BLUE, label="Culex sightings")
ax.plot(prof.index, prof["new_pools"], "-s", ms=3, color=ORANGE, label="WNV-positive pools (DSHS report week)")
ax.set_xlabel("MMWR week"); ax.set_ylabel("share of 2019-2025 total")
ax.set_title("Texas seasonal timing, 2019-2025 pooled", fontsize=10); ax.legend(frameon=False, fontsize=8)
fig.savefig(os.path.join(FIG, "tx_seasonal.png")); plt.close(fig)


# ----------------------------------------------------------------- fixed-effects machinery
def demean(v, groups, iters=30):
    """residualize v on several sets of dummies by alternating projections"""
    v = v.astype(float).copy()
    for _ in range(iters):
        prev = v.copy()
        for gr in groups:
            v -= pd.Series(v).groupby(gr).transform("mean").values
        if np.abs(v - prev).max() < 1e-9:
            break
    return v


def fe_beta(df, y, x, fe_cols, cluster="county"):
    """OLS slope of y on x net of the fixed effects in fe_cols, county-clustered SE"""
    df = df.dropna(subset=[y, x])
    codes = [pd.factorize(df[c])[0] for c in fe_cols]
    yt, xt = demean(df[y].values, codes), demean(df[x].values, codes)
    sxx = (xt ** 2).sum()
    if sxx == 0:
        return None
    b = (xt * yt).sum() / sxx
    e = yt - b * xt
    gi = pd.factorize(df[cluster])[0]
    s = pd.Series(xt * e).groupby(gi).sum().values
    G, n = len(s), len(df)
    k = sum(len(np.unique(c)) for c in codes)
    se = np.sqrt((s ** 2).sum() * G / (G - 1) * (n - 1) / max(n - k, 1)) / sxx
    p = 2 * (1 - stats.t.cdf(abs(b / se), G - 1))
    return {"beta": float(b), "se": float(se), "p": float(p), "n": int(n), "counties": int(G)}


def rnd(r):
    return None if r is None else {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}


# ----------------------------------------------------------------- 1. county-week lag test
panel_all["yearweek"] = panel_all["year"] * 100 + panel_all["week"]
panel_all["county_year"] = panel_all["county"] + panel_all["year"].astype(str)


def lagged(df, x, k):
    """signal from k weeks earlier, same county"""
    return df.groupby("county")[x].shift(k)


def season_rows(df):
    return df[df["week"].between(*SEASON)]


lag_tab = {}
for ykey, ylab in [("y_pools", "pools"), ("y_human", "human")]:
    for x in ["x_share", "x_obsr", "x_culex"]:
        for k in LAGS:
            df = panel_all.copy()
            df["xl"] = lagged(df, x, k)
            df = season_rows(df)
            r = fe_beta(df, ykey, "xl", ["county", "yearweek"])
            lag_tab[f"{ylab}|{x}|lag{k}"] = rnd(r)
R["lag_county_week"] = lag_tab

# variant: within-county-season (county x year FE + week FE): only timing inside a county's season identifies
alt = {}
for k in LAGS:
    df = panel_all.copy()
    df["xl"] = lagged(df, "x_share", k)
    alt[f"lag{k}"] = rnd(fe_beta(season_rows(df), "y_pools", "xl", ["county_year", "week"]))
R["lag_within_county_season_share"] = alt

# anomaly specification: county-specific seasonal shape (county x week-of-year) plus year-week effects, so
# only a county's year-to-year deviation from its own usual week identifies the effect
panel_all["county_week"] = panel_all["county"] + "_" + panel_all["week"].astype(str)
anom = {}
for ykey, ylab in [("y_pools", "pools"), ("y_human", "human")]:
    for x in ["x_share", "x_obsr"]:
        for k in LAGS:
            df = panel_all.copy()
            df["xl"] = lagged(df, x, k)
            anom[f"{ylab}|{x}|lag{k}"] = rnd(fe_beta(season_rows(df), ykey, "xl", ["county_week", "yearweek"]))
R["lag_anomaly_spec"] = anom


def placebo(k, fe_cols, n=N_PERM):
    """pair each county's signal with a different year's pools; one-sided p = share of shuffles at least as large"""
    base = panel_all.copy()
    base["xl"] = lagged(base, "x_share", k)
    base = season_rows(base).dropna(subset=["xl", "y_pools"])
    obs_b = fe_beta(base, "y_pools", "xl", fe_cols)["beta"]
    years = sorted(base["year"].unique())
    xp_all = base.set_index(["county", "week", "year"])["xl"]
    perm = []
    for _ in range(n):
        while True:
            pi = dict(zip(years, rng.permutation(years)))
            if all(pi[y] != y for y in years):
                break
        key = pd.MultiIndex.from_arrays([base["county"], base["week"], base["year"].map(pi)])
        r = fe_beta(base.assign(xp=xp_all.reindex(key).values), "y_pools", "xp", fe_cols)
        if r:
            perm.append(r["beta"])
    perm = np.array(perm)
    return {"beta_observed": round(obs_b, 4), "perm_mean": round(float(perm.mean()), 4), "perm_sd": round(float(perm.std()), 4),
            "p_one_sided": round(float((np.sum(perm >= obs_b) + 1) / (len(perm) + 1)), 4),
            "z_vs_perm": round(float((obs_b - perm.mean()) / perm.std()), 2), "n_perm": int(len(perm))}


R["placebo"] = {"county_yearweek_lag2": placebo(2, ["county", "yearweek"]),
                "county_yearweek_lag5": placebo(5, ["county", "yearweek"]),
                "anomaly_spec_lag2": placebo(2, ["county_week", "yearweek"]),
                "anomaly_spec_lag5": placebo(5, ["county_week", "yearweek"])}

# ----------------------------------------------------------------- 2. county-year season test
early = (14, 26)
ey = panel_all[panel_all["week"].between(*early)].groupby(["county", "year"])[["mosq_obs", "all_obs", "mosq_observers", "culex_obs"]].sum()
ey["x_share"] = share(ey["mosq_obs"], ey["all_obs"])
ey["x_obsr"] = share(ey["mosq_observers"], ey["all_obs"])
ey["x_culex"] = share(ey["culex_obs"], ey["all_obs"])
sy = season_rows(panel_all).groupby(["county", "year"])[["new_pools", "new_human_total"]].sum()
cy = ey.join(sy).reset_index()
cy["y_pools"] = np.log1p(cy["new_pools"])
cy["y_human"] = np.log1p(cy["new_human_total"])
season = {}
for ykey in ["y_pools", "y_human"]:
    for x in ["x_share", "x_obsr", "x_culex"]:
        season[f"{ykey}|{x}"] = rnd(fe_beta(cy, ykey, x, ["county", "year"]))
R["season_county_year"] = season
R["season_county_year"]["n_county_years"] = int(len(cy))

# year-level statewide: early mosquito share vs the year's total pools
st = panel_all[panel_all["week"].between(*early)].groupby("year")[["mosq_obs", "all_obs"]].sum()
st["share_per_1000"] = 1000 * st["mosq_obs"] / st["all_obs"]
st = st.join(panel_all.groupby("year")["new_pools"].sum())
rho, p = stats.spearmanr(st["share_per_1000"], st["new_pools"])
R["year_level"] = {"spearman_rho": round(float(rho), 3), "p": round(float(p), 4), "n": int(len(st)),
                   "table": st.round(3).reset_index().to_dict("records")}

# leave-one-year-out forecast of each county's season pools
def loyo(x):
    eb, em = [], []
    for y0 in sorted(cy["year"].unique()):
        tr, te = cy[cy["year"] != y0], cy[cy["year"] == y0]
        my, mx = tr.groupby("county")["y_pools"].mean(), tr.groupby("county")[x].mean()
        td = tr.assign(yd=tr["y_pools"] - tr["county"].map(my), xd=tr[x] - tr["county"].map(mx))
        b = np.polyfit(td["xd"], td["yd"], 1)[0]
        te = te[te["county"].isin(my.index)]
        base_ = te["county"].map(my)
        pred = base_ + b * (te[x] - te["county"].map(mx))
        eb += list((te["y_pools"] - base_).abs()); em += list((te["y_pools"] - pred).abs())
    b_, m_ = float(np.mean(eb)), float(np.mean(em))
    return {"mae_baseline": round(b_, 4), "mae_with_signal": round(m_, 4),
            "improvement_pct": round(100 * (b_ - m_) / b_, 2), "n": len(eb)}


R["out_of_sample"] = {x: loyo(x) for x in ["x_share", "x_obsr", "x_culex"]}

# ----------------------------------------------------------------- 3. robustness of the anomaly-spec result
ROB_K = 4                  # centre of the 3-6 week plateau; chosen after seeing the lag curve


def spec_beta(df, k=ROB_K, x="x_share", y="y_pools", fe_cols=("county_week", "yearweek")):
    df = df.copy()
    df["yearweek"] = df["year"] * 100 + df["week"]
    df["county_week"] = df["county"] + "_" + df["week"].astype(str)
    df["xl"] = lagged(df, x, k)
    return rnd(fe_beta(season_rows(df), y, "xl", list(fe_cols)))


top_pools = list(tot.loc[list(tested)].sort_values("pools", ascending=False).index)
rob = {"baseline": spec_beta(panel_all)}
for name, minobs in [("min_obs_100k", 100000), ("min_obs_5k", 5000)]:
    cs = sorted(set(tot.index[tot["all_obs"] >= minobs]) & tested)
    rob[name] = {"counties": len(cs), **spec_beta(t[t["county"].isin(cs)])}
rob["drop_top4_pool_counties"] = {"dropped": top_pools[:4], **spec_beta(panel_all[~panel_all["county"].isin(top_pools[:4])])}
rob["drop_top10_pool_counties"] = {"dropped": top_pools[:10], **spec_beta(panel_all[~panel_all["county"].isin(top_pools[:10])])}
rob["only_top10_pool_counties"] = spec_beta(panel_all[panel_all["county"].isin(top_pools[:10])])
for y0 in sorted(panel_all["year"].unique()):
    rob[f"drop_year_{y0}"] = spec_beta(panel_all[panel_all["year"] != y0])
rob["county_year_fe_variant"] = spec_beta(panel_all, fe_cols=("county_year", "week"))
rob["culex_signal"] = spec_beta(panel_all, x="x_culex")
R["robustness_anomaly_lag4"] = rob

# incremental value: does the sighting signal still matter once the county's own pool history is known?
def fe_multi(df, y, xs, fe_cols):
    """slope on xs[0] with the other xs as controls, net of fixed effects, county-clustered SE"""
    df = df.dropna(subset=[y] + xs)
    codes = [pd.factorize(df[c])[0] for c in fe_cols]
    yt = demean(df[y].values, codes)
    X = np.column_stack([demean(df[c].values, codes) for c in xs])
    b, *_ = np.linalg.lstsq(X, yt, rcond=None)
    e = yt - X @ b
    bread = np.linalg.inv(X.T @ X)
    gi = pd.factorize(df["county"])[0]
    meat = sum(np.outer(s, s) for s in
               (pd.DataFrame(X * e[:, None]).groupby(gi).sum().values))
    G, n = len(np.unique(gi)), len(df)
    kk = sum(len(np.unique(c)) for c in codes) + len(xs)
    V = bread @ meat @ bread * G / (G - 1) * (n - 1) / max(n - kk, 1)
    se = float(np.sqrt(V[0, 0]))
    return {"beta": round(float(b[0]), 4), "se": round(se, 4),
            "p": round(float(2 * (1 - stats.t.cdf(abs(b[0] / se), G - 1))), 4), "n": int(n),
            "beta_cum_pools_control": round(float(b[1]), 4) if len(xs) > 1 else None}


pc = panel_all.copy()
pc["yearweek"] = pc["year"] * 100 + pc["week"]
pc["county_week"] = pc["county"] + "_" + pc["week"].astype(str)
pc["cum_before"] = np.log1p(pc.groupby(["county", "year"])["cum_pools"].shift(1))     # pools reported through last week
ctrl = {}
for k in range(1, 9):
    pc["xl"] = lagged(pc, "x_share", k)
    pc["cum_k"] = pc.groupby("county")["cum_before"].shift(k - 1)                       # known when the signal was read
    s = season_rows(pc)
    ctrl[f"lag{k}"] = {"without_control": fe_multi(s, "y_pools", ["xl"], ["county_week", "yearweek"]),
                       "with_cum_pools_control": fe_multi(s, "y_pools", ["xl", "cum_k"], ["county_week", "yearweek"])}
R["lag_with_pool_history_control"] = ctrl

# effect size: one SD of the (windowed) mosquito-share signal, in plain terms
tmp = season_rows(panel_all).dropna(subset=["x_share"])
sd_x = float(tmp["x_share"].std())
b4 = anom[f"pools|x_share|lag{ROB_K}"]["beta"]
R["effect_size"] = {"sd_x_share": round(sd_x, 3), "beta_lag4": b4, "log1p_pools_change_per_sd": round(b4 * sd_x, 4),
                    "mean_new_pools_per_county_week": round(float(tmp["new_pools"].mean()), 3),
                    "pct_county_weeks_with_new_pool": round(100 * float((tmp["new_pools"] > 0).mean()), 1)}

# ----------------------------------------------------------------- lag figure
fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.2), sharey=True)
for ax, (x, title) in zip(axes, [("x_share", "Mosquito share of observations"), ("x_obsr", "Unique mosquito observers per observation")]):
    ks = list(LAGS)
    for key, col, lab, off in [("pools", BLUE, "positive mosquito pools", -0.08), ("human", ORANGE, "human cases", 0.08)]:
        b = np.array([anom[f"{key}|{x}|lag{k}"]["beta"] for k in ks])
        se = np.array([anom[f"{key}|{x}|lag{k}"]["se"] for k in ks])
        ax.errorbar(np.array(ks) + off, b, yerr=1.96 * se, fmt="o", ms=4, color=col, capsize=2, label=lab)
    ax.axhline(0, color=GREY, lw=.6)
    ax.set_xlabel("weeks the signal leads the outcome"); ax.set_title(title, fontsize=9)
axes[0].legend(frameon=False, fontsize=8)
axes[0].set_ylabel("effect on log(1+new count)\n(county-week + year-week effects, 95% CI)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "tx_lags.png")); plt.close(fig)

json.dump(R, open(os.path.join(D, "tx_results.json"), "w"), indent=1, default=float)
print(json.dumps(R["effect_size"], indent=1))
print("\nCONTROLLING FOR THE COUNTY'S OWN POOL HISTORY (anomaly spec)")
for k, v in R["lag_with_pool_history_control"].items():
    print(" ", k, "no control:", v["without_control"]["beta"], v["without_control"]["p"],
          "| with control:", v["with_cum_pools_control"]["beta"], v["with_cum_pools_control"]["se"], v["with_cum_pools_control"]["p"],
          "| pool-history coef:", v["with_cum_pools_control"]["beta_cum_pools_control"])
print("\nROBUSTNESS (anomaly spec, lag 4, pools ~ share)")
for k, v in R["robustness_anomaly_lag4"].items():
    print(" ", k, {a: b for a, b in v.items() if a in ("beta", "se", "p", "n", "counties")})
print(json.dumps(R["placebo"], indent=1))
for name, tab_ in [("county + yearweek FE, pools, share", lambda k: lag_tab[f"pools|x_share|lag{k}"]),
                   ("county x year + week FE, pools, share", lambda k: alt[f"lag{k}"]),
                   ("ANOMALY county x week + yearweek FE, pools, share", lambda k: anom[f"pools|x_share|lag{k}"]),
                   ("ANOMALY, pools, observers", lambda k: anom[f"pools|x_obsr|lag{k}"]),
                   ("ANOMALY, human, share", lambda k: anom[f"human|x_share|lag{k}"])]:
    print("\n" + name)
    for k in LAGS:
        r = tab_(k); print(" ", k, r["beta"], r["se"], r["p"], r["n"])
