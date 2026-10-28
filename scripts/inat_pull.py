"""
iNaturalist pull for the Mosquito Sightings / West Nile project.
Works for California (CA) or Texas (TX).

What it does
  1. Looks up taxon IDs for Culicidae (mosquitoes) and Culex.
  2. Finds iNaturalist's place for the state and each county, using a point
     inside each county (from Census boundaries) so spelling doesn't matter.
  3. Downloads every Culicidae observation in the state for the chosen years
     (all quality grades, so you can filter later).
  4. Downloads daily counts of ALL iNaturalist observations per county
     (the baseline for normal activity).
  5. Builds a county x CDC (MMWR) week table. For CA, pass the CDPH case file
     with --wnv to merge human cases.

Outputs (in ./data_<state> by default, e.g. ./data_tx)
  inat_mosquito_observations_<st>.csv    one row per mosquito observation
  inat_all_obs_daily_by_county_<st>.csv  one row per county-day
  inat_county_week_<st>.csv              the analysis table
  cache/                                 saved progress; reruns pick up where they stopped

Usage
  pip install requests pandas
  python inat_pull.py --state TX --start 2019 --end 2025
  python inat_pull.py --state CA --start 2015 --end 2025 --wnv wnv_human_cases.csv

Runtime: Texas has 254 counties, so expect about 1 to 1.5 hours the first time
(mostly the baseline step, one request per second). California takes 15 to 30 minutes.
Follows iNaturalist's guidance: ~1 request/second, 200 results per page,
id_above paging, backoff on HTTP 429.
"""

import argparse
import json
import os
import time
from datetime import date, timedelta

import pandas as pd
import requests

VERSION = "5 (multi-state)"
API = "https://api.inaturalist.org/v1"
HEADERS = {"User-Agent": "MosquitoSightingsResearch/0.2 (Cornell student project)"}
PAUSE = 1.1

STATES = {
    "CA": {"name": "California"},
    "TX": {"name": "Texas"},
}

# A point inside each county, by state (from Census TIGER via us-atlas).
COUNTY_POINTS = json.loads(r'''{"CA":{"San Benito":[36.5825,-121.0221],"Shasta":[40.7331,-122.0179],"Sierra":[39.5991,-120.5099],"Fresno":[36.7569,-119.6306],"Amador":[38.4678,-120.7257],"El Dorado":[38.8006,-120.5211],"Monterey":[36.3558,-121.4916],"San Bernardino":[34.8507,-116.1193],"Tehama":[40.1339,-122.2146],"Nevada":[39.2529,-121.0247],"Tuolumne":[38.0361,-119.8844],"Santa Barbara":[34.73,-120.0308],"Trinity":[40.6571,-123.1382],"Stanislaus":[37.558,-120.9772],"Colusa":[39.1477,-122.208],"Alpine":[38.6353,-119.8441],"Siskiyou":[41.4925,-122.5415],"Del Norte":[41.6918,-123.9086],"Inyo":[36.6156,-117.5009],"San Joaquin":[37.8969,-121.2457],"Imperial":[33.0265,-115.354],"Calaveras":[38.131,-120.6677],"Butte":[39.7181,-121.588],"Madera":[37.2612,-119.6441],"Ventura":[34.6273,-119.1228],"Orange":[33.6811,-117.7371],"Kern":[35.3073,-118.7198],"Contra Costa":[37.9192,-121.9762],"Kings":[36.1599,-119.7617],"Mendocino":[39.3844,-123.4432],"Sonoma":[38.4819,-122.857],"Yuba":[39.2774,-121.4371],"Merced":[37.2011,-120.6069],"Tulare":[36.2389,-118.8199],"Sacramento":[38.3947,-121.2763],"Solano":[38.2925,-121.897],"San Francisco":[37.7479,-122.4405],"Mariposa":[37.5275,-119.9524],"Marin":[38.052,-122.7404],"Sutter":[39.0141,-121.7105],"Lassen":[40.4457,-120.5288],"Riverside":[33.7339,-116.0183],"Plumas":[40.051,-120.7999],"Santa Clara":[37.1873,-121.718],"Placer":[39.0188,-120.7379],"San Luis Obispo":[35.3073,-120.3435],"San Mateo":[37.3961,-122.3086],"Mono":[38.0906,-119.0341],"San Diego":[33.0437,-116.7002],"Napa":[38.4914,-122.3226],"Alameda":[37.6976,-121.8818],"Yolo":[38.6228,-121.8542],"Lake":[39.1477,-122.7856],"Los Angeles":[34.2657,-118.1445],"Glenn":[39.6125,-122.44],"Santa Cruz":[37.0769,-122.0376],"Humboldt":[40.7138,-123.9036],"Modoc":[41.6029,-120.7244]},"TX":{"Walker":[30.744,-95.5909],"Brooks":[26.919,-98.2034],"Deaf Smith":[35.0689,-102.6056],"Hockley":[33.6072,-102.3427],"Shackelford":[32.7357,-99.3537],"Tom Green":[31.3058,-100.4019],"Wise":[33.307,-97.6538],"Cottle":[34.0299,-100.2836],"Taylor":[32.2984,-99.8895],"Armstrong":[34.9649,-101.3575],"Hutchinson":[35.8403,-101.3548],"Irion":[31.3053,-100.9806],"Mitchell":[32.3057,-100.9205],"Yoakum":[33.1098,-102.8282],"Childress":[34.437,-100.2062],"Gillespie":[30.3926,-98.9476],"Dickens":[33.6158,-100.7787],"La Salle":[28.4468,-99.0989],"Hansford":[36.2784,-101.3539],"Potter":[35.4018,-101.8941],"Falls":[31.2504,-96.9377],"Hamilton":[31.7628,-98.0692],"Runnels":[31.8878,-99.9768],"Hale":[34.0715,-101.8259],"San Jacinto":[30.6323,-95.1853],"Duval":[27.708,-98.5175],"Moore":[35.8416,-101.8932],"Martin":[32.3053,-101.9515],"Limestone":[31.4677,-96.5199],"Zavala":[28.8669,-99.762],"Cochran":[33.6974,-102.8304],"Fisher":[32.743,-100.4018],"Baylor":[33.5664,-99.214],"Panola":[32.1777,-94.3108],"Jack":[33.2017,-98.172],"Coryell":[31.3921,-97.7998],"Lipscomb":[36.278,-100.2726],"Donley":[34.9649,-100.8146],"Atascosa":[28.9863,-98.5455],"Briscoe":[34.5298,-101.2085],"Goliad":[28.6483,-97.466],"Gray":[35.4013,-100.8137],"Hall":[34.5307,-100.6809],"Brown":[31.8895,-98.9968],"Gaines":[32.7408,-102.6353],"Bell":[31.0524,-97.5134],"Sutton":[30.4983,-100.5382],"San Augustine":[31.3822,-94.1795],"King":[33.6166,-100.2565],"Uvalde":[29.357,-99.762],"Wilson":[29.168,-98.0815],"Schleicher":[30.8964,-100.5382],"Swisher":[34.5303,-101.7352],"Motley":[34.0741,-100.7796],"Lynn":[33.1759,-101.816],"Menard":[30.8256,-99.7997],"Nolan":[32.3044,-100.4054],"Sterling":[31.8959,-101.0442],"Garza":[33.1828,-101.2973],"Callahan":[32.298,-99.3735],"Stonewall":[33.1802,-100.2547],"Franklin":[33.2009,-95.2171],"Van Zandt":[32.5836,-95.8362],"Hunt":[33.1094,-96.0794],"Madison":[30.9527,-95.9388],"Hood":[32.4153,-97.8215],"Mason":[30.8161,-99.2249],"Gregg":[32.4788,-94.8438],"Hemphill":[35.7501,-100.2708],"Castro":[34.5303,-102.2629],"Upton":[31.4248,-102.0448],"Jim Hogg":[27.0985,-98.7113],"Lavaca":[29.3596,-96.9008],"Midland":[31.8689,-102.0314],"Kimble":[30.3952,-99.7091],"Carson":[35.4009,-101.3548],"Palo Pinto":[32.7559,-98.3176],"Lamb":[34.0689,-102.3517],"Eastland":[32.2799,-98.8146],"Kerr":[29.9914,-99.3053],"Lubbock":[33.6098,-101.8205],"Erath":[32.1605,-98.1926],"Upshur":[32.7301,-94.9273],"Blanco":[30.2478,-98.403],"Andrews":[32.3044,-102.6379],"Jeff Davis":[30.716,-104.1659],"Knox":[33.5664,-99.734],"Llano":[30.6903,-98.6861],"Wheeler":[35.303,-100.2695],"Crosby":[33.6141,-101.3],"McMullen":[28.4073,-98.5677],"Roberts":[35.8403,-100.8146],"Kent":[33.1841,-100.7778],"Comal":[29.8007,-98.2969],"Zapata":[26.9448,-99.168],"Aransas":[28.2355,-96.9524],"Orange":[30.0563,-93.8728],"Winkler":[31.8264,-103.0633],"Jones":[32.74,-99.8787],"Sherman":[36.2776,-101.8932],"Chambers":[29.7071,-94.5252],"Hardin":[30.3123,-94.3865],"Hill":[31.9952,-97.1201],"Williamson":[30.6817,-97.6341],"Kinney":[29.3789,-100.4081],"Calhoun":[28.3953,-96.5981],"Ellis":[32.2946,-96.7699],"Clay":[33.6514,-98.2007],"Kaufman":[32.6794,-96.2984],"Foard":[33.9457,-99.7624],"Collin":[33.1644,-96.568],"Comanche":[31.9608,-98.5922],"Concho":[31.2774,-99.8583],"Crane":[31.4905,-102.5409],"Jasper":[30.7379,-94.005],"Gonzales":[29.4459,-97.5046],"Jefferson":[29.8763,-94.1051],"Dimmit":[28.4468,-99.7557],"Bowie":[33.4354,-94.3946],"Brewster":[29.823,-103.0814],"Dallam":[36.1865,-102.6021],"Cass":[33.1304,-94.3473],"Fayette":[29.8772,-96.929],"Burnet":[30.7225,-98.1941],"Refugio":[28.2952,-97.2014],"Ochiltree":[36.2784,-100.8155],"Matagorda":[28.7883,-95.951],"Throckmorton":[33.1772,-99.2129],"Anderson":[31.8096,-95.6915],"Caldwell":[29.851,-97.6267],"Medina":[29.3866,-99.1087],"Glasscock":[31.8689,-101.5208],"Kenedy":[26.8984,-97.7206],"Wood":[32.8117,-95.3989],"Liberty":[30.2564,-94.8906],"Travis":[30.2818,-97.7594],"Montgomery":[30.3286,-95.5116],"Howard":[32.3057,-101.4355],"Bastrop":[30.1375,-97.3098],"Trinity":[31.1026,-95.1264],"McLennan":[31.5549,-97.2037],"McCulloch":[31.2486,-99.3476],"Brazoria":[29.2161,-95.4666],"Bailey":[34.1244,-102.83],"Houston":[31.2267,-95.4625],"San Saba":[31.2143,-98.8288],"Presidio":[29.9497,-104.2391],"Nacogdoches":[31.5339,-94.5905],"Washington":[30.2161,-96.4323],"Cooke":[33.5754,-97.2158],"Wilbarger":[34.142,-99.2138],"Leon":[31.2946,-96.0091],"Karnes":[28.9446,-97.9017],"Randall":[34.9654,-101.8968],"DeWitt":[29.0555,-97.3663],"Ward":[31.4823,-103.1315],"Wichita":[33.9491,-98.6879],"Willacy":[26.4645,-97.6863],"Shelby":[31.7868,-94.1305],"Kendall":[29.8595,-98.7135],"Morris":[33.1257,-94.7348],"Hopkins":[33.1094,-95.5859],"Rockwall":[32.9118,-96.4069],"Hudspeth":[31.3036,-105.4084],"Smith":[32.4449,-95.2673],"Hidalgo":[26.3907,-98.1898],"Brazos":[30.6499,-96.3811],"Burleson":[30.5068,-96.6262],"Delta":[33.3533,-95.7263],"Crockett":[30.6903,-101.3772],"Lee":[30.2779,-96.8856],"Robertson":[31.0335,-96.5165],"Mills":[31.4677,-98.6428],"Loving":[31.8204,-103.569],"Sabine":[31.3809,-93.8514],"Austin":[29.8565,-96.2842],"Reagan":[31.3041,-101.5235],"Reeves":[31.3784,-103.6431],"Tarrant":[32.7726,-97.2917],"Hays":[30.1972,-98.093],"San Patricio":[28.0151,-97.5093],"Oldham":[35.4018,-102.6039],"Pecos":[30.7255,-102.6691],"Real":[29.8724,-99.8096],"Dallas":[32.6811,-96.7782],"Nueces":[27.7798,-97.6623],"Red River":[33.5505,-95.0282],"Coleman":[31.8277,-99.4584],"Rusk":[32.109,-94.7896],"Terrell":[30.168,-102.1628],"Tyler":[30.7929,-94.3854],"Victoria":[28.7965,-96.9371],"Webb":[27.7368,-99.2938],"Wharton":[29.2909,-96.2228],"Kleberg":[27.419,-97.7143],"Harris":[29.8445,-95.3951],"Hardeman":[34.343,-99.736],"Newton":[30.7109,-93.7506],"Navarro":[32.0433,-96.467],"Montague":[33.6587,-97.7323],"Frio":[28.8691,-99.107],"Harrison":[32.5368,-94.3718],"Lamar":[33.6617,-95.5832],"Floyd":[34.0732,-101.3027],"Starr":[26.5011,-98.786],"Stephens":[32.7348,-98.8351],"Parmer":[34.5303,-102.7851],"Marion":[32.7713,-94.3503],"Val Verde":[29.7685,-101.0801],"El Paso":[31.6584,-106.1606],"Coke":[31.8891,-100.5301],"Cameron":[26.1171,-97.5735],"Waller":[29.9867,-96.0006],"Bee":[28.4232,-97.7109],"Titus":[33.2198,-94.9697],"Camp":[32.9981,-94.9985],"Grimes":[30.5451,-95.9905],"Scurry":[32.7456,-100.916],"Jim Wells":[27.7063,-98.0858],"Ector":[31.8689,-102.5428],"Rains":[32.8568,-95.8058],"Parker":[32.7765,-97.8041],"Young":[33.2017,-98.6889],"Polk":[30.8032,-94.8882],"Angelina":[31.2787,-94.6334],"Live Oak":[28.3979,-98.135],"Dawson":[32.7421,-101.9479],"Haskell":[33.179,-99.7315],"Terry":[33.1738,-102.3347],"Archer":[33.6506,-98.6878],"Henderson":[32.1915,-95.8181],"Bosque":[31.9101,-97.6657],"Hartley":[35.8974,-102.6021],"Collingsworth":[34.8893,-100.2702],"Somervell":[32.1455,-97.8243],"Colorado":[29.5971,-96.5271],"Freestone":[31.7267,-96.1636],"Borden":[32.7443,-101.4319],"Culberson":[31.4595,-104.4926],"Galveston":[29.3995,-95.0388],"Lampasas":[31.3204,-98.248],"Edwards":[29.9609,-100.348],"Guadalupe":[29.622,-97.9485],"Cherokee":[31.7744,-95.174],"Fannin":[33.5771,-96.1191],"Fort Bend":[29.5262,-95.7433],"Milam":[30.7916,-96.9704],"Bexar":[29.3441,-98.5109],"Johnson":[32.4342,-97.3509],"Denton":[33.2945,-97.113],"Bandera":[29.7509,-99.2254],"Maverick":[28.6427,-100.2856],"Jackson":[28.9923,-96.5816],"Grayson":[33.5711,-96.6635]}}''')


# ---------------------------------------------------------------- helpers

def get(path, params=None, tries=6):
    wait = 5
    for _ in range(tries):
        time.sleep(PAUSE)
        r = requests.get(f"{API}/{path}", params=params, headers=HEADERS, timeout=60)
        if r.status_code == 200:
            return r.json()
        if r.status_code in (429, 500, 502, 503, 504):
            print(f"  HTTP {r.status_code}, waiting {wait}s")
            time.sleep(wait)
            wait = min(wait * 2, 300)
            continue
        r.raise_for_status()
    raise RuntimeError(f"Gave up on {path} {params}")


def mmwr_week(d):
    def week1_start(y):
        jan1 = date(y, 1, 1)
        k = (jan1.weekday() + 1) % 7
        s = jan1 - timedelta(days=k)
        return s + timedelta(days=7) if k > 3 else s
    y = d.year
    if d >= week1_start(y + 1):
        y += 1
    elif d < week1_start(y):
        y -= 1
    return y, (d - week1_start(y)).days // 7 + 1


def load_json(path, default):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return default


def save_json(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f)


# ---------------------------------------------------------------- taxa

def find_taxon(name, rank):
    for t in get("taxa", {"q": name, "rank": rank, "per_page": 30})["results"]:
        if t["name"].lower() == name.lower() and t["rank"] == rank:
            print(f"  {name}: taxon_id {t['id']}")
            return t["id"]
    raise RuntimeError(f"Could not find taxon {name} ({rank})")


# ---------------------------------------------------------------- places

def find_state(state_name):
    for q in (state_name, f"{state_name}, US"):
        for p in get("places/autocomplete", {"q": q, "per_page": 50})["results"]:
            if p.get("admin_level") == 10 and (p.get("name") or "").lower() == state_name.lower():
                return p["id"], p.get("display_name")
        for p in get("places/autocomplete", {"q": q, "per_page": 50})["results"]:
            if (p.get("display_name") or "").lower() == f"{state_name}, us".lower():
                return p["id"], p.get("display_name")
    raise RuntimeError(f"Could not find the place for {state_name}")


def find_county_by_point(county, lat, lng, state_id):
    ids = set()
    for radius in (3, 10, 30):
        for o in get("observations", {"lat": lat, "lng": lng, "radius": radius,
                                      "per_page": 10})["results"]:
            ids.update(o.get("place_ids") or [])
        if ids:
            break
    if not ids:
        return None, None
    ids = sorted(ids)
    places = []
    for i in range(0, len(ids), 50):
        places += get("places/" + ",".join(map(str, ids[i:i + 50])))["results"]
    names = {county.lower(), f"{county} county".lower()}
    level20 = [p for p in places if p.get("admin_level") == 20
               and state_id in (p.get("ancestor_place_ids") or [state_id])]
    for p in level20:
        if (p.get("name") or "").strip().lower() in names:
            return p["id"], p.get("display_name")
    if level20:
        print(f"  {county}: county places near the point were "
              f"{[(p['id'], p.get('name')) for p in level20]}, none matched by name")
    return None, None


def find_county_by_name(county, abbr, state_id):
    target = f"{county} County, {abbr}, US".lower()
    names = {county.lower(), f"{county} county".lower()}
    seen = {}
    for q in (f"{county} County", county):
        for p in get("places/autocomplete", {"q": q, "per_page": 50})["results"]:
            seen[p["id"]] = p
            if (p.get("display_name") or "").strip().lower() == target:
                return p["id"], p["display_name"]
    cands = [p for p in seen.values()
             if state_id in (p.get("ancestor_place_ids") or [])
             and (p.get("name") or "").strip().lower() in names]
    if cands:
        cands.sort(key=lambda p: (p.get("admin_level") == 20, p.get("bbox_area") or 0), reverse=True)
        return cands[0]["id"], cands[0].get("display_name")
    print(f"  No name match for {county}. Candidates:")
    for p in list(seen.values())[:6]:
        print(f"    id={p['id']} name={p.get('name')!r} display={p.get('display_name')!r} "
              f"admin_level={p.get('admin_level')}")
    return None, None


def get_places(abbr, cache_dir):
    path = os.path.join(cache_dir, "places.json")
    places = load_json(path, {})
    state_name = STATES[abbr]["name"]
    if not places.get("_state"):
        pid, disp = find_state(state_name)
        places["_state"] = pid
        print(f"  {state_name}: place_id {pid} ({disp})")
        save_json(path, places)
    for county, (lat, lng) in sorted(COUNTY_POINTS[abbr].items()):
        if places.get(county):
            continue
        pid, disp = find_county_by_point(county, lat, lng, places["_state"])
        if not pid:
            pid, disp = find_county_by_name(county, abbr, places["_state"])
        places[county] = pid
        print(f"  {county}: {'place_id ' + str(pid) + ' (' + str(disp) + ')' if pid else 'NOT FOUND'}")
        save_json(path, places)
    missing = [c for c in COUNTY_POINTS[abbr] if not places.get(c)]
    if missing:
        print(f"  WARNING: no place found for {missing}. You can add IDs by hand in {path}.")
    return places


# ---------------------------------------------------------------- mosquito records

def pull_mosquitoes(taxon_id, state_place, county_by_place, culex_id, start, end, out_csv, cache_dir):
    state_path = os.path.join(cache_dir, "mosq_state.json")
    state = load_json(state_path, {"id_above": 0})
    rows = []
    if os.path.exists(out_csv) and state["id_above"] > 0:
        rows = pd.read_csv(out_csv).to_dict("records")
        print(f"  Resuming after id {state['id_above']} ({len(rows)} rows so far)")
    params = {"taxon_id": taxon_id, "place_id": state_place,
              "d1": f"{start}-01-01", "d2": f"{end}-12-31",
              "per_page": 200, "order_by": "id", "order": "asc"}
    total = get("observations", {**params, "per_page": 0})["total_results"]
    print(f"  {total} Culicidae observations, {start}-{end}")
    while True:
        data = get("observations", {**params, "id_above": state["id_above"]})["results"]
        if not data:
            break
        for o in data:
            taxon = o.get("taxon") or {}
            counties = [county_by_place[p] for p in (o.get("place_ids") or []) if p in county_by_place]
            lat = lon = None
            if o.get("location"):
                lat, lon = (float(x) for x in o["location"].split(","))
            rows.append({
                "id": o["id"], "observed_on": o.get("observed_on"),
                "created_at": (o.get("created_at") or "")[:10],
                "user_id": (o.get("user") or {}).get("id"),
                "quality_grade": o.get("quality_grade"),
                "taxon_id": taxon.get("id"), "taxon_name": taxon.get("name"),
                "taxon_rank": taxon.get("rank"),
                "is_culex": culex_id in (taxon.get("ancestor_ids") or []) or taxon.get("id") == culex_id,
                "county": counties[0] if counties else None,
                "latitude": lat, "longitude": lon,
                "positional_accuracy": o.get("positional_accuracy"),
                "obscured": o.get("obscured"), "geoprivacy": o.get("geoprivacy"),
            })
        state["id_above"] = data[-1]["id"]
        pd.DataFrame(rows).to_csv(out_csv, index=False)
        save_json(state_path, state)
        print(f"  {len(rows)} rows downloaded")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- baseline

def pull_baseline(abbr, places, start, end, out_csv, cache_dir):
    path = os.path.join(cache_dir, "baseline.json")
    done = load_json(path, {})
    counties = sorted(COUNTY_POINTS[abbr])
    for i, c in enumerate(counties, 1):
        pid = places.get(c)
        if not pid:
            continue
        for y in range(start, end + 1):
            key = f"{c}|{y}"
            if key in done:
                continue
            res = get("observations/histogram", {"place_id": pid, "d1": f"{y}-01-01",
                                                 "d2": f"{y}-12-31", "date_field": "observed",
                                                 "interval": "day"})
            done[key] = res["results"].get("day", {})
            save_json(path, done)
        print(f"  [{i}/{len(counties)}] {c} done")
    rows = [{"county": k.split("|")[0], "date": d, "all_obs": n}
            for k, days in done.items() for d, n in days.items()]
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    return df


# ---------------------------------------------------------------- county-week table

def add_mmwr(df, col):
    d = pd.to_datetime(df[col], errors="coerce")
    df = df[d.notna()].copy()
    yw = [mmwr_week(x.date()) for x in d[d.notna()]]
    df["Year"] = [a for a, _ in yw]
    df["Week"] = [b for _, b in yw]
    return df


def build_table(abbr, mosq, base, start, end, wnv_path, out_csv):
    m = add_mmwr(mosq.dropna(subset=["county"]), "observed_on")
    m_agg = m.groupby(["county", "Year", "Week"]).agg(
        mosq_obs=("id", "count"), culex_obs=("is_culex", "sum"),
        mosq_observers=("user_id", "nunique"),
        mosq_research_grade=("quality_grade", lambda s: (s == "research").sum()),
    ).reset_index()
    b = add_mmwr(base, "date").groupby(["county", "Year", "Week"]).all_obs.sum().reset_index()
    weeks = sorted({mmwr_week(date(start, 1, 1) + timedelta(days=i))
                    for i in range((date(end, 12, 31) - date(start, 1, 1)).days + 1)})
    grid = pd.DataFrame([(c, y, w) for c in sorted(COUNTY_POINTS[abbr]) for (y, w) in weeks
                         if start <= y <= end], columns=["county", "Year", "Week"])
    t = grid.merge(b, how="left").merge(m_agg, how="left")
    cols = ["all_obs", "mosq_obs", "culex_obs", "mosq_observers", "mosq_research_grade"]
    t[cols] = t[cols].fillna(0).astype(int)
    t["mosq_share"] = (t.mosq_obs / t.all_obs).where(t.all_obs > 0)
    if wnv_path:
        w = pd.read_csv(wnv_path)
        w.columns = ["Year", "Week", "county", "wnv_cases"]
        t = t.merge(w, how="left", on=["county", "Year", "Week"])
        t["wnv_cases"] = t.wnv_cases.fillna(0).astype(int)
    t.to_csv(out_csv, index=False)
    return t


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", choices=sorted(STATES), default="CA")
    ap.add_argument("--start", type=int, default=2019)
    ap.add_argument("--end", type=int, default=2025)
    ap.add_argument("--wnv", help="CA only: path to CDPH wnv_human_cases.csv")
    ap.add_argument("--out", help="output folder (default data_<state>)")
    args = ap.parse_args()

    abbr = args.state
    st = abbr.lower()
    out = args.out or f"data_{st}"
    cache = os.path.join(out, "cache")
    os.makedirs(cache, exist_ok=True)
    print(f"inat_pull.py version {VERSION}, state {abbr}, {args.start}-{args.end}")

    print("Step 1: taxa")
    culicidae = find_taxon("Culicidae", "family")
    culex = find_taxon("Culex", "genus")

    print("Step 2: places")
    places = get_places(abbr, cache)
    county_by_place = {pid: c for c, pid in places.items() if pid and c != "_state"}

    print("Step 3: mosquito observations")
    mosq = pull_mosquitoes(culicidae, places["_state"], county_by_place, culex, args.start, args.end,
                           os.path.join(out, f"inat_mosquito_observations_{st}.csv"), cache)

    print("Step 4: all-observation counts by county and day")
    base = pull_baseline(abbr, places, args.start, args.end,
                         os.path.join(out, f"inat_all_obs_daily_by_county_{st}.csv"), cache)

    print("Step 5: county-week table")
    t = build_table(abbr, mosq, base, args.start, args.end, args.wnv,
                    os.path.join(out, f"inat_county_week_{st}.csv"))

    print("\nSummary")
    print(f"  counties with a place: {len(county_by_place)} of {len(COUNTY_POINTS[abbr])}")
    print(f"  mosquito observations: {len(mosq)}  (no county: {mosq.county.isna().sum()})")
    print(f"  county-weeks with >=1 mosquito obs: {(t.mosq_obs > 0).mean():.1%}")
    print(f"\nDone. Files are in ./{out}/")


if __name__ == "__main__":
    main()
