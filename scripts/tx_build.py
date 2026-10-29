"""
Merge the Texas iNaturalist county-week table with the parsed DSHS weekly WNV counts.

Reads   data_tx/inat_county_week_tx.csv     (scripts/inat_pull.py --state TX)
        data_tx/tx_wnv_county_week.csv      (scripts/tx_parse_reports.py)
Writes  data_tx/tx_panel.csv                one row per county x MMWR week, 2019-2025

DSHS weekly counts are new reports that week (difference of year-to-date totals), so they date
from when DSHS reported a pool, not when it was collected. A week's outcome is set to missing when
that week or the one before has no DSHS report (the difference would be lumped), and counties not
listed in a report had no activity, so missing counties are zeros.
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data_tx")

FIX = {"Hutchison": "Hutchinson", "Kleburg": "Kleberg"}      # spelling errors in the DSHS tables

inat = pd.read_csv(os.path.join(D, "inat_county_week_tx.csv")).rename(columns={"Year": "year", "Week": "week"})
w = pd.read_csv(os.path.join(D, "tx_wnv_county_week.csv"))
w["county"] = w["county"].replace(FIX)
w = w[w["year"].between(2019, 2025)]

unknown = sorted(set(w["county"]) - set(inat["county"]))
assert not unknown, f"DSHS counties not in the iNaturalist table: {unknown}"

# weeks with a DSHS report, and weeks whose difference is clean (previous week also reported)
have = {(y, k) for y, k in zip(w["year"], w["week"])}
# week-1 reports start the year from zero; any other first report of a year holds a lump of earlier weeks
ok = {(y, k) for (y, k) in have if (y, k - 1) in have or k == 1}

cols = ["new_pools", "new_avian", "new_equine", "new_chicken", "new_wnf", "new_wnnd", "new_human_total", "cum_pools"]
w = w.groupby(["county", "year", "week"], as_index=False)[cols].sum()
t = inat.merge(w, how="left", on=["county", "year", "week"])
t["reported"] = [(y, k) in have for y, k in zip(t["year"], t["week"])]
t["clean"] = [(y, k) in ok for y, k in zip(t["year"], t["week"])]
t[cols] = t[cols].fillna(0)
# year-to-date totals carry forward through weeks where a county is unlisted
t = t.sort_values(["county", "year", "week"])
t["cum_pools"] = t.groupby(["county", "year"])["cum_pools"].cummax()
for c in cols[:-1]:
    t.loc[~t["clean"], c] = np.nan
t.to_csv(os.path.join(D, "tx_panel.csv"), index=False)

print(f"{len(t)} rows, {t['county'].nunique()} counties, {t['year'].min()}-{t['year'].max()}")
print(f"clean outcome weeks: {int(t['clean'].mean() * 100)}% of rows")
print(t.groupby("year")[["new_pools", "new_human_total", "mosq_obs", "all_obs"]].sum().round(0).to_string())
