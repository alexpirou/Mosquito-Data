"""Turn raw pulls into analysis tables.

Outputs (data/processed/):
  mosquito_obs.csv      one row per verifiable US mosquito observation, with state + county FIPS
  state_month.csv       state x month: mosquito obs, observers, Culex obs, insect / crane-fly / all-taxa obs
  national_month.csv    national month series joined to CDC cases by onset month
  state_year.csv        state x year: early-season (Apr-Jul) and full-year crowd signals + CDC cases
  county.csv            county, pooled 2016-2025: crowd signals + CDC cases + population
"""
import json, os
import numpy as np
import pandas as pd
from shapely.geometry import Polygon, MultiPolygon, Point
from shapely.strtree import STRtree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, INAT, CDC = (os.path.join(ROOT, "data", "raw"), os.path.join(ROOT, "data", "raw", "inat"),
                  os.path.join(ROOT, "data", "raw", "cdc"))
OUT = os.path.join(ROOT, "data", "processed"); os.makedirs(OUT, exist_ok=True)
YEARS = list(range(2016, 2026))          # CDC data finalized through 2025
EARLY = [4, 5, 6, 7]                     # Apr-Jul: before the Aug-Sep peak in human onsets

STATES = {  # name: (fips, abbr)
    "Alabama": ("01", "AL"), "Alaska": ("02", "AK"), "Arizona": ("04", "AZ"), "Arkansas": ("05", "AR"),
    "California": ("06", "CA"), "Colorado": ("08", "CO"), "Connecticut": ("09", "CT"), "Delaware": ("10", "DE"),
    "District of Columbia": ("11", "DC"), "Florida": ("12", "FL"), "Georgia": ("13", "GA"), "Hawaii": ("15", "HI"),
    "Idaho": ("16", "ID"), "Illinois": ("17", "IL"), "Indiana": ("18", "IN"), "Iowa": ("19", "IA"),
    "Kansas": ("20", "KS"), "Kentucky": ("21", "KY"), "Louisiana": ("22", "LA"), "Maine": ("23", "ME"),
    "Maryland": ("24", "MD"), "Massachusetts": ("25", "MA"), "Michigan": ("26", "MI"), "Minnesota": ("27", "MN"),
    "Mississippi": ("28", "MS"), "Missouri": ("29", "MO"), "Montana": ("30", "MT"), "Nebraska": ("31", "NE"),
    "Nevada": ("32", "NV"), "New Hampshire": ("33", "NH"), "New Jersey": ("34", "NJ"), "New Mexico": ("35", "NM"),
    "New York": ("36", "NY"), "North Carolina": ("37", "NC"), "North Dakota": ("38", "ND"), "Ohio": ("39", "OH"),
    "Oklahoma": ("40", "OK"), "Oregon": ("41", "OR"), "Pennsylvania": ("42", "PA"), "Rhode Island": ("44", "RI"),
    "South Carolina": ("45", "SC"), "South Dakota": ("46", "SD"), "Tennessee": ("47", "TN"), "Texas": ("48", "TX"),
    "Utah": ("49", "UT"), "Vermont": ("50", "VT"), "Virginia": ("51", "VA"), "Washington": ("53", "WA"),
    "West Virginia": ("54", "WV"), "Wisconsin": ("55", "WI"), "Wyoming": ("56", "WY"),
}
FIPS2ABBR = {f: a for f, a in STATES.values()}


# ------------------------------------------------------------------ county polygons (TopoJSON)
def load_counties():
    topo = json.load(open(os.path.join(RAW, "us_counties_2021_topo.json")))
    tr = topo.get("transform")
    arcs = []
    for arc in topo["arcs"]:
        if tr:
            x = y = 0; pts = []
            for dx, dy in arc:
                x += dx; y += dy
                pts.append((x * tr["scale"][0] + tr["translate"][0], y * tr["scale"][1] + tr["translate"][1]))
        else:
            pts = [tuple(p) for p in arc]
        arcs.append(pts)

    def ring(idx):
        pts = []
        for i in idx:
            a = arcs[i] if i >= 0 else arcs[~i][::-1]
            pts.extend(a if not pts else a[1:])
        return pts

    obj = topo["objects"]["counties"]
    geoms, fips = [], []
    for g in obj["geometries"]:
        props = g.get("properties", {})
        gid = props.get("GEOID") or props.get("geoid") or g.get("id")
        if g["type"] == "Polygon":
            rings = [ring(r) for r in g["arcs"]]
            shp = Polygon(rings[0], rings[1:])
        elif g["type"] == "MultiPolygon":
            shp = MultiPolygon([Polygon(ring(p[0]), [ring(r) for r in p[1:]]) for p in g["arcs"]])
        else:
            continue
        if not shp.is_valid:
            shp = shp.buffer(0)
        geoms.append(shp); fips.append(str(gid).zfill(5))
    return geoms, fips


def assign_fips(df):
    geoms, fips = load_counties()
    tree = STRtree(geoms)
    out = []
    for lat, lng in zip(df["lat"], df["lng"]):
        if pd.isna(lat):
            out.append(None); continue
        p = Point(lng, lat)
        hits = tree.query(p, predicate="intersects")
        out.append(fips[hits[0]] if len(hits) else None)
    return out


# ------------------------------------------------------------------ iNaturalist mosquito observations
def build_obs():
    rows = [json.loads(l) for l in open(os.path.join(INAT, "culicidae_us.jsonl"))]
    df = pd.DataFrame(rows)
    places = json.load(open(os.path.join(INAT, "places.json")))
    state_place = {int(k): STATES[v["name"]][0] for k, v in places.items()
                   if v["admin_level"] == 10 and v["name"] in STATES and 1 in (v["ancestor_place_ids"] or [])}
    county_place = {int(k) for k, v in places.items() if v["admin_level"] == 20}

    df["poly_fips"] = assign_fips(df)
    df["place_state"] = [next((state_place[p] for p in pids if p in state_place), None) for pids in df["place_ids"]]
    df["county_place"] = [next((p for p in pids if p in county_place), None) for pids in df["place_ids"]]

    # iNat county place -> FIPS by majority vote over un-obscured observations
    ok = df[(~df["obscured"].fillna(False)) & df["poly_fips"].notna() & df["county_place"].notna()]
    cmap = ok.groupby("county_place")["poly_fips"].agg(lambda s: s.value_counts().index[0]).to_dict()
    # county: exact point for un-obscured records; for obscured ones prefer iNat's own county membership
    df["fips"] = np.where(df["obscured"].fillna(False) & df["county_place"].map(cmap).notna(),
                          df["county_place"].map(cmap), df["poly_fips"])
    df["state_fips"] = df["fips"].str[:2].where(df["fips"].notna(), df["place_state"])
    df = df[df["state_fips"].isin(FIPS2ABBR)]
    df = df[~df["captive"].fillna(False)]
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df[df["date"].notna()]
    df["year"], df["month"] = df["date"].dt.year, df["date"].dt.month
    df["state"] = df["state_fips"].map(FIPS2ABBR)
    keep = ["id", "date", "year", "month", "state", "state_fips", "fips", "user", "taxon", "rank", "culex",
            "qg", "obscured", "lat", "lng"]
    df[keep].to_csv(os.path.join(OUT, "mosquito_obs.csv"), index=False)
    json.dump({str(k): v for k, v in cmap.items()}, open(os.path.join(OUT, "county_place_to_fips.json"), "w"))
    print("obs:", len(df), "rows;", df["fips"].notna().mean().round(3), "with county")
    return df, cmap


# ------------------------------------------------------------------ CDC
def load_cdc():
    cy = pd.read_csv(os.path.join(CDC, "county_year.csv"), dtype={"County": str})
    cy.columns = ["year", "fips", "activity", "cases", "neuro", "blood_donor", "notes"]
    cy = cy[cy["year"].astype(str).str.isdigit()].copy()
    cy["year"] = cy["year"].astype(int)
    cy["fips"] = cy["fips"].str.zfill(5)
    for c in ["cases", "neuro", "blood_donor"]:
        cy[c] = pd.to_numeric(cy[c], errors="coerce").fillna(0)
    cy["state_fips"] = cy["fips"].str[:2]
    nm = pd.read_csv(os.path.join(CDC, "national_onset_month.csv"))
    nm.columns = ["type", "year", "month", "cases"]
    nm = nm[nm["year"].astype(str).str.isdigit()].copy()
    nm["year"] = nm["year"].astype(int)
    nm["month"] = pd.to_datetime(nm["month"], format="%b").dt.month
    nm = nm.pivot_table(index=["year", "month"], columns="type", values="cases").reset_index()
    nm.columns = ["year", "month", "wnv_cases", "wnv_neuro"]
    return cy, nm


# ------------------------------------------------------------------ tables
def hist_frame():
    h = json.load(open(os.path.join(INAT, "state_month_hist.json")))
    rows = []
    for pid, d in h.items():
        if d["name"] not in STATES:
            continue
        sf = STATES[d["name"]][0]
        for series in ["insecta", "culicidae", "tipulidae", "culex", "all"]:
            for ts, n in d.get(series, {}).items():
                t = pd.Timestamp(ts)
                rows.append((sf, t.year, t.month, series, n))
    df = pd.DataFrame(rows, columns=["state_fips", "year", "month", "series", "n"])
    return df.pivot_table(index=["state_fips", "year", "month"], columns="series", values="n",
                          fill_value=0).reset_index()


def main():
    obs, cmap = build_obs()
    cy, nm = load_cdc()
    hist = hist_frame().rename(columns={"insecta": "insect_obs", "tipulidae": "cranefly_obs",
                                        "all": "all_obs", "culicidae": "mosq_obs_api", "culex": "culex_obs_api"})

    # ---------- state x month
    g = obs.groupby(["state_fips", "year", "month"])
    sm = pd.DataFrame({"mosq_obs": g.size(), "mosq_observers": g["user"].nunique(),
                       "culex_obs": g["culex"].sum(),
                       "mosq_rg": g["qg"].apply(lambda s: (s == "research").sum())}).reset_index()
    sm = hist.merge(sm, on=["state_fips", "year", "month"], how="left").fillna(0)
    sm["state"] = sm["state_fips"].map(FIPS2ABBR)
    sm.to_csv(os.path.join(OUT, "state_month.csv"), index=False)

    # ---------- national x month
    nat_users = obs.groupby(["year", "month"])["user"].nunique().rename("mosq_observers")
    nat = sm.groupby(["year", "month"])[["mosq_obs", "culex_obs", "insect_obs", "cranefly_obs", "all_obs",
                                         "mosq_rg"]].sum()
    nat = nat.join(nat_users).reset_index().merge(nm, on=["year", "month"], how="left")
    nat.to_csv(os.path.join(OUT, "national_month.csv"), index=False)

    # ---------- state x year
    def window(df, months, tag):
        d = df[df["month"].isin(months)]
        a = d.groupby(["state_fips", "year"])[["mosq_obs", "culex_obs", "insect_obs", "cranefly_obs",
                                               "all_obs", "mosq_rg"]].sum()
        u = obs[obs["month"].isin(months)].groupby(["state_fips", "year"])["user"].nunique().rename("mosq_observers")
        a = a.join(u).fillna(0)
        return a.add_suffix("_" + tag)

    sy = window(sm, EARLY, "early").join(window(sm, range(1, 13), "full")).reset_index()
    sy = sy[sy["year"].isin(YEARS)]
    cases = cy.groupby(["state_fips", "year"])[["cases", "neuro"]].sum().reset_index()
    sy = sy.merge(cases, on=["state_fips", "year"], how="left").fillna({"cases": 0, "neuro": 0})
    sy["state"] = sy["state_fips"].map(FIPS2ABBR)
    sy.to_csv(os.path.join(OUT, "state_year.csv"), index=False)

    # ---------- county, pooled 2016-2025
    o = obs[obs["year"].isin(YEARS) & obs["fips"].notna()]
    cg = o.groupby("fips")
    county = pd.DataFrame({"mosq_obs": cg.size(), "mosq_observers": cg["user"].nunique(),
                           "culex_obs": cg["culex"].sum()})
    ins = json.load(open(os.path.join(INAT, "county_year_insecta.json")))
    ins_f = {}
    for pid, d in ins.items():
        f = cmap.get(int(pid))
        if f is None:
            continue
        ins_f[f] = ins_f.get(f, 0) + sum(n for ts, n in d.items() if int(ts[:4]) in YEARS)
    county = county.join(pd.Series(ins_f, name="insect_obs"), how="outer").fillna(0)
    c = cy[cy["year"].isin(YEARS)]
    county = county.join(c.groupby("fips")[["cases", "neuro"]].sum(), how="left").fillna(0)
    act = c[c["activity"].notna() & (c["activity"] != "Data unavailable")]
    county = county.join(act.groupby("fips")["year"].nunique().rename("years_any_activity"),
                         how="left").fillna({"years_any_activity": 0})
    pop = pd.read_csv(os.path.join(RAW, "county_pop_plotly.csv"), dtype={"FIPS": str})
    pop["FIPS"] = pop["FIPS"].str.zfill(5)
    county = county.join(pop.set_index("FIPS")["TOT_POP"].rename("pop"), how="left")
    county.index.name = "fips"
    county = county.reset_index()
    county["state"] = county["fips"].str[:2].map(FIPS2ABBR)
    county.to_csv(os.path.join(OUT, "county.csv"), index=False)
    print("tables written:", len(sm), "state-months;", len(sy), "state-years;", len(county), "counties")


if __name__ == "__main__":
    main()
