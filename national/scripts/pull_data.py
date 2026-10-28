"""Pull raw data for the iNaturalist mosquito / West Nile virus project.

Steps (each is resumable; finished outputs are skipped on re-run):
  cdc       CDC ArboNET historic WNV files (county-year, state-year, national onset month)
  topo      Census county boundaries (TopoJSON served by cdc.gov) for point-in-county
  taxa      iNaturalist taxon ids used in the study
  mosq      every verifiable US mosquito (Culicidae) observation, 2015-2025
  places    lookup of the state/county "standard places" referenced by those observations
  hist      monthly histograms per US state for insects, mosquitoes, crane flies, all taxa
  county    yearly insect counts for every county that has at least one mosquito record

Usage:  python pull_data.py [step ...]   (no args = all steps in order)
"""
import json, os, subprocess, sys, threading, time
import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
INAT, CDC = os.path.join(RAW, "inat"), os.path.join(RAW, "cdc")
os.makedirs(INAT, exist_ok=True); os.makedirs(CDC, exist_ok=True)

API = "https://api.inaturalist.org"
D1, D2 = "2015-01-01", "2025-12-31"
UA = "CS5306-INFO5306 class project (Cornell) - WNV hidden signals"
PAUSE = 1.0  # per worker; with at most 4 workers we stay near iNaturalist's ~60 requests/minute guidance
_local = threading.local()


def session():
    if not hasattr(_local, "s"):
        _local.s = requests.Session(); _local.s.headers["User-Agent"] = UA
    return _local.s


def pmap(fn, items, workers=3):
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(workers) as ex:
        return list(ex.map(fn, items))


def get_json(url, params=None, tries=6):
    for i in range(tries):
        try:
            r = session().get(url, params=params, timeout=90)
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(15 * (i + 1)); continue
            r.raise_for_status()
            time.sleep(PAUSE)
            return r.json()
        except requests.RequestException as e:
            print("  retry", i, e, flush=True); time.sleep(10 * (i + 1))
    raise RuntimeError(f"failed: {url} {params}")


# ---------------------------------------------------------------- CDC
CDC_BASE = "https://www.cdc.gov/wcms/vizdata/live/ncezid_dvbd/WNV/"
CDC_FILES = {
    "county_year.csv": "wnv_hist_hum_nonhum_yearly.csv",
    "county_cumulative.csv": "wnv_hist_hum_nonhum_cumulative.csv",
    "county_incidence.csv": "wnv_hum_historic_County%20Inc.csv",
    "state_year.csv": "wnv_hum_historic%20-%20Historic%20Sta.csv",
    "national_onset_month.csv": "wnv_hum_historic_Year%20and%20Month.csv",
    "national_year.csv": "wnv_hum_historic%20-%20Yearly%20Data.csv",
    "national_summary.csv": "wnv_hum_historic_Data%20Bites.csv",
}


def curl(url, path):
    """cdc.gov serves these files to plain curl; python-requests gets a 403"""
    subprocess.run(["curl", "-sS", "--fail", "-m", "120", "-o", path, url], check=True)


def step_cdc():
    for local, remote in CDC_FILES.items():
        path = os.path.join(CDC, local)
        if os.path.exists(path):
            continue
        curl(CDC_BASE + remote, path)
        print("cdc:", local, os.path.getsize(path), "bytes")


def step_topo():
    path = os.path.join(RAW, "us_counties_2021_topo.json")
    if os.path.exists(path):
        return
    url = "https://www.cdc.gov/TemplatePackage/contrib/data/county-topography/cb_2021_us_county_20m.json"
    curl(url, path)
    print("topo:", os.path.getsize(path), "bytes")


# ---------------------------------------------------------------- taxa
def step_taxa():
    path = os.path.join(INAT, "taxa.json")
    if os.path.exists(path):
        return
    want = {"Culicidae": "family", "Culex": "genus", "Tipulidae": "family", "Insecta": "class"}
    out = {}
    for name, rank in want.items():
        res = get_json(f"{API}/v1/taxa", {"q": name, "rank": rank})["results"]
        hit = [t for t in res if t["name"] == name][0]
        out[name] = hit["id"]
    json.dump(out, open(path, "w"), indent=1)
    print("taxa:", out)


def taxa():
    return json.load(open(os.path.join(INAT, "taxa.json")))


# ---------------------------------------------------------------- mosquito observations
FIELDS = ("(id:!t,observed_on:!t,user:(id:!t),taxon:(id:!t,rank:!t,ancestor_ids:!t),"
          "place_ids:!t,quality_grade:!t,location:!t,obscured:!t,captive:!t)")


# Four date slices pulled concurrently (~50 requests/min in total, under iNaturalist's limit)
SLICES = [("2015-01-01", "2021-12-31"), ("2022-01-01", "2023-12-31"),
          ("2024-01-01", "2024-12-31"), ("2025-01-01", "2025-12-31")]


def step_mosq():
    from concurrent.futures import ThreadPoolExecutor
    final = os.path.join(INAT, "culicidae_us.jsonl")
    if os.path.exists(os.path.join(INAT, "culicidae_us.done")):
        return
    with ThreadPoolExecutor(len(SLICES)) as ex:
        list(ex.map(lambda s: pull_slice(*s), SLICES))
    seen = set()
    with open(final, "w") as out:
        for d1, d2 in SLICES:
            for line in open(os.path.join(INAT, f"culicidae_{d1[:4]}_{d2[:4]}.jsonl")):
                i = json.loads(line)["id"]
                if i not in seen:
                    seen.add(i); out.write(line)
    open(os.path.join(INAT, "culicidae_us.done"), "w").write(str(len(seen)))
    print("mosq: merged", len(seen), "observations", flush=True)


def pull_slice(d1, d2):
    tag = f"{d1[:4]}_{d2[:4]}"
    path = os.path.join(INAT, f"culicidae_{tag}.jsonl")
    done = path + ".done"
    if os.path.exists(done):
        return
    tx = taxa()
    last = 0
    if os.path.exists(path):
        for line in open(path):
            last = max(last, json.loads(line)["id"])
    n_total, n = None, 0
    with open(path, "a") as f:
        while True:
            url = (f"{API}/v2/observations?taxon_id={tx['Culicidae']}&place_id=1&verifiable=true"
                   f"&d1={d1}&d2={d2}&order_by=id&order=asc&per_page=200&id_above={last}&fields={FIELDS}")
            js = get_json(url)
            n_total = n_total or js["total_results"]
            res = js["results"]
            n += len(res)
            if not res:
                break
            for o in res:
                lat, lng = (o.get("location") or ",").split(",") if o.get("location") else (None, None)
                anc = (o.get("taxon") or {}).get("ancestor_ids") or []
                f.write(json.dumps({
                    "id": o["id"], "date": o.get("observed_on"), "user": (o.get("user") or {}).get("id"),
                    "taxon": (o.get("taxon") or {}).get("id"), "rank": (o.get("taxon") or {}).get("rank"),
                    "culex": tx["Culex"] in anc, "qg": o.get("quality_grade"),
                    "lat": float(lat) if lat else None, "lng": float(lng) if lng else None,
                    "obscured": o.get("obscured"), "captive": o.get("captive"),
                    "place_ids": o.get("place_ids") or [],
                }) + "\n")
            last = res[-1]["id"]
            f.flush()
            print(f"mosq {tag}: {n}/{n_total}", flush=True)
    open(done, "w").write("ok")


# ---------------------------------------------------------------- places
def step_places():
    path = os.path.join(INAT, "places.json")
    if os.path.exists(path):
        return
    ids = set()
    for line in open(os.path.join(INAT, "culicidae_us.jsonl")):
        ids.update(json.loads(line)["place_ids"])
    ids = sorted(ids)
    print("places: looking up", len(ids), "place ids", flush=True)
    def lookup(i):
        js = get_json(f"{API}/v1/places/{','.join(map(str, ids[i:i + 100]))}")
        return [{k: p.get(k) for k in ("id", "name", "display_name", "admin_level", "ancestor_place_ids")}
                for p in js["results"] if p.get("admin_level") in (10, 20)]
    keep = {p["id"]: p for chunk in pmap(lookup, range(0, len(ids), 100)) for p in chunk}
    json.dump(keep, open(path, "w"))
    print("places: kept", len(keep))


def us_places():
    pl = json.load(open(os.path.join(INAT, "places.json")))
    us = [p for p in pl.values() if p["ancestor_place_ids"] and 1 in p["ancestor_place_ids"]]
    states = [p for p in us if p["admin_level"] == 10]
    counties = [p for p in us if p["admin_level"] == 20]
    return states, counties


# ---------------------------------------------------------------- state places (by name)
STATE_NAMES = ["Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut", "Delaware",
               "District of Columbia", "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa",
               "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota",
               "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire", "New Jersey",
               "New Mexico", "New York", "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon",
               "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah",
               "Vermont", "Virginia", "Washington", "West Virginia", "Wisconsin", "Wyoming"]


def step_states():
    path = os.path.join(INAT, "state_places.json")
    if os.path.exists(path):
        return
    out = []
    for name in STATE_NAMES:
        res = get_json(f"{API}/v1/places/autocomplete", {"q": name, "per_page": 20})["results"]
        hit = [p for p in res if p["name"] == name and p.get("admin_level") == 10
               and 1 in (p.get("ancestor_place_ids") or [])]
        if not hit:
            print("states: no match for", name, [(p["name"], p.get("admin_level")) for p in res[:5]]); continue
        out.append({k: hit[0].get(k) for k in ("id", "name", "display_name", "admin_level", "ancestor_place_ids")})
    json.dump(out, open(path, "w"))
    print("states:", len(out), "matched")


# ---------------------------------------------------------------- state-month histograms
def step_hist():
    path = os.path.join(INAT, "state_month_hist.json")
    out = json.load(open(path)) if os.path.exists(path) else {}
    tx = taxa()
    states = json.load(open(os.path.join(INAT, "state_places.json")))
    series = {"insecta": tx["Insecta"], "culicidae": tx["Culicidae"], "tipulidae": tx["Tipulidae"],
              "culex": tx["Culex"], "all": None}
    for st in states:
        key = str(st["id"])
        out.setdefault(key, {"name": st["name"]})
        for s, tid in series.items():
            if s in out[key]:
                continue
            params = {"place_id": st["id"], "interval": "month", "date_field": "observed",
                      "d1": D1, "d2": D2, "verifiable": "true"}
            if tid:
                params["taxon_id"] = tid
            out[key][s] = get_json(f"{API}/v1/observations/histogram", params)["results"]["month"]
        json.dump(out, open(path, "w"))
        print("hist:", st["name"], flush=True)


# ---------------------------------------------------------------- county-year insect denominators
def step_county():
    path = os.path.join(INAT, "county_year_insecta.json")
    out = json.load(open(path)) if os.path.exists(path) else {}
    tx = taxa()
    _, counties = us_places()
    todo = [c for c in counties if str(c["id"]) not in out]
    print("county:", len(todo), "counties to query", flush=True)

    def one(c):
        params = {"place_id": c["id"], "taxon_id": tx["Insecta"], "interval": "year",
                  "date_field": "observed", "d1": D1, "d2": D2, "verifiable": "true"}
        return str(c["id"]), get_json(f"{API}/v1/observations/histogram", params)["results"]["year"]
    for k in range(0, len(todo), 150):        # checkpoint every 150 counties
        out.update(dict(pmap(one, todo[k:k + 150])))
        json.dump(out, open(path, "w"))
        print(f"county: {len(out)}/{len(counties)}", flush=True)
    json.dump(out, open(path, "w"))
    print("county: done", len(out))


# ---------------------------------------------------------------- genus composition (context only)
def step_genera():
    path = os.path.join(INAT, "genus_counts.json")
    if os.path.exists(path):
        return
    out = {}
    for g in ["Aedes", "Culex", "Anopheles", "Psorophora", "Culiseta", "Toxorhynchites", "Coquillettidia", "Uranotaenia"]:
        res = get_json(f"{API}/v1/taxa", {"q": g, "rank": "genus"})["results"]
        hit = [t for t in res if t["name"] == g and "Culicidae" in str(t.get("ancestry", "")) or t["name"] == g]
        tid = hit[0]["id"]
        n = get_json(f"{API}/v1/observations", {"taxon_id": tid, "place_id": 1, "verifiable": "true",
                                                "d1": "2016-01-01", "d2": D2, "per_page": 0})["total_results"]
        out[g] = {"taxon_id": tid, "obs_2016_2025": n}
    json.dump(out, open(path, "w"), indent=1)
    print("genera:", out)



STEPS = {"cdc": step_cdc, "topo": step_topo, "taxa": step_taxa, "states": step_states, "mosq": step_mosq,
         "places": step_places, "hist": step_hist, "county": step_county, "genera": step_genera}

if __name__ == "__main__":
    for name in (sys.argv[1:] or list(STEPS)):
        print("==", name, flush=True)
        STEPS[name]()
