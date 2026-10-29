"""
Parse the Texas DSHS weekly arbovirus reports (PDFs from tx_download_reports.py)
into a county-week table of West Nile activity.

Reads   tx_reports_manifest.csv and the PDFs in tx_reports/
Writes  data_tx/tx_wnv_county_cumulative.csv   one row per report x county (year-to-date, as printed)
        data_tx/tx_wnv_county_week.csv         weekly new counts (difference of consecutive reports)
        data_tx/tx_parse_log.txt               checks and warnings

Notes
  - The report year and week come from the first-page header ("Week #35 (ending August 31, 2019)"),
    not the file name: some file names carry the wrong year (e.g. 01Jan2022Week52 is a 2021 report).
  - Table 3's grid is ragged and its columns move around between years, so the PDF words are read
    with their x positions. The last four non-"H" letter headers (M, A, E/V, SC) belong to the WNV
    block; WNF, WNND, PVD and TOTAL follow it. Each number goes to the nearest of those 8 anchors.
  - Counties with no activity are not listed, so missing = 0. Early-season reports with no Table 3
    rows are valid empty weeks.
  - Where two reports exist for one (year, week), the one whose file name says "Corrected" wins,
    otherwise the later week-ending date.
"""
import os
import re
from datetime import datetime

import pandas as pd
import pdfplumber

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data_tx")
os.makedirs(OUT, exist_ok=True)
LOG = []


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


HEAD_RE = re.compile(r"Week\s*#?\s*(\d{1,2})\s*\(ending\s+([A-Za-z]+)\s+(\d{1,2}),\s*(\d{4})\)", re.I)
FIELDS = ["pools", "avian", "equine", "chicken", "wnf", "wnnd", "pvd", "human_total"]
LETTERS = {"M", "A", "E", "V", "SC"}


def clean(w):
    return re.sub(r"[^A-Za-z0-9]", "", w)


def xc(w):
    return (w["x0"] + w["x1"]) / 2


def lines_of(words, tol=3):
    rows = []
    for w in sorted(words, key=lambda w: (w["top"], w["x0"])):
        if rows and abs(rows[-1][0] - w["top"]) <= tol:
            rows[-1][1].append(w)
        else:
            rows.append([w["top"], [w]])
    return [sorted(r[1], key=lambda w: w["x0"]) for r in rows]


def find_anchors(page_words):
    """x centres of the 8 WNV columns, from the header words of the page's Table 3"""
    ls = lines_of(page_words)
    for i, ln in enumerate(ls):
        toks = [clean(w["text"]) for w in ln]
        if "WNF" in toks and "WNND" in toks:
            sub = ln
            # the letter row sits just above the WNF row
            for j in range(i - 1, max(i - 4, -1), -1):
                lt = [(clean(w["text"]).upper().rstrip("0123456789"), w) for w in ls[j]]
                lt = [(t, w) for t, w in lt if t in LETTERS]
                if len(lt) >= 4:
                    wnv = [w for t, w in lt if t != "H"][-4:]
                    hdr = {clean(w["text"]).upper(): w for w in sub}
                    tot = [w for w in sub if clean(w["text"]).upper() == "TOTAL"]
                    pvd = [w for w in sub if clean(w["text"]).upper().startswith("PVD")]
                    if len(wnv) == 4 and pvd and tot and "WNF" in hdr and "WNND" in hdr:
                        # TOTAL: the one right of PVD (2019 style also has COUNTY TOTAL further right)
                        t = sorted([w for w in tot if w["x0"] > pvd[0]["x1"] - 1], key=lambda w: w["x0"])[0]
                        return [xc(w) for w in wnv] + [xc(hdr["WNF"]), xc(hdr["WNND"]), xc(pvd[0]), xc(t)], ln[0]["top"]
    return None, None


def parse_pdf(path):
    with pdfplumber.open(path) as pdf:
        text1 = pdf.pages[0].extract_text() or ""
        m = HEAD_RE.search(text1)
        if not m:
            return None, "no header"
        week, mon, day, year = int(m.group(1)), m.group(2), int(m.group(3)), int(m.group(4))
        end = datetime.strptime(f"{mon[:3]} {day} {year}", "%b %d %Y").date()
        # MMWR year = year in the report title; the week-ending date can fall in early January of the next year
        ty = re.match(r"\s*(\d{4})\s+DSHS", text1)
        if ty and abs(int(ty.group(1)) - year) <= 1:
            year = int(ty.group(1))
        mm = re.search(r"^West Nile\s+(\d+)", text1, re.M)
        t1 = int(mm.group(1)) if mm else None

        rows, anchors, in_t3, seen_t3 = [], None, False, False
        for page in pdf.pages:
            t = page.extract_text() or ""
            if re.search(r"Table\s*3\.", t):
                in_t3 = seen_t3 = True
            elif re.search(r"Figure\s*1\.|Table\s*[124]\.", t):
                in_t3 = False
            if not in_t3:
                continue
            words = page.extract_words(keep_blank_chars=False, use_text_flow=False)
            a, htop = find_anchors(words)
            if a:
                anchors = a
            if anchors is None:
                continue
            gap = min(b - c for c, b in zip(anchors, anchors[1:])) / 2
            body = lines_of([w for w in words if (htop is None or w["top"] > htop + 1)])
            # some early-season reports print the county name on the line below its numbers
            merged, skip = [], False
            for i, ln in enumerate(body):
                if skip:
                    skip = False
                    continue
                nxt = body[i + 1] if i + 1 < len(body) else None
                if (nxt and re.fullmatch(r"\d+", ln[0]["text"]) and not any(re.search(r"\d", w["text"]) for w in nxt)
                        and re.fullmatch(r"[A-Z][A-Za-z .'\-]*", " ".join(w["text"] for w in nxt))
                        and not " ".join(w["text"] for w in nxt).lower().startswith(("total", "of ", "m ", "number", "report"))
                        and 0 < nxt[0]["top"] - ln[0]["top"] < 15):
                    ln = nxt + ln
                    skip = True
                merged.append(ln)
            for ln in merged:
                texts = [w["text"] for w in ln]
                first = texts[0]
                if first.lower().startswith(("total", "of", "m", "reports")) and not re.fullmatch(r"[A-Z][a-z]+", first):
                    if first.lower().startswith("total"):
                        break
                    continue
                num_i = next((i for i, w in enumerate(ln) if re.fullmatch(r"\d+", w["text"])), None)
                if num_i is None or num_i == 0:
                    continue
                name = " ".join(texts[:num_i])
                if re.fullmatch(r"[A-Za-z][A-Za-z .'\-]*", name) is None or len(name) > 20:
                    continue
                if name.lower().startswith(("number", "total", "report", "of reports")):
                    continue
                rec = dict.fromkeys(FIELDS, 0)
                for w in ln[num_i:]:
                    if not re.fullmatch(r"\d+", w["text"]):
                        continue
                    d = [abs(xc(w) - c) for c in anchors]
                    k = d.index(min(d))
                    if d[k] <= gap:
                        rec[FIELDS[k]] += int(w["text"])
                rows.append({"county": name, **rec})
    return {"year": year, "week": week, "week_end": end, "t1_pools": t1,
            "df": pd.DataFrame(rows, columns=["county"] + FIELDS), "seen_t3": seen_t3}, None


def main():
    man = pd.read_csv(os.path.join(ROOT, "tx_reports_manifest.csv"))
    reps, bad = [], []
    for n, f in enumerate(man["file"], 1):
        path = os.path.join(ROOT, f)
        try:
            r, err = parse_pdf(path)
        except Exception as e:
            r, err = None, f"{type(e).__name__}: {e}"
        if r is None:
            bad.append((f, err))
        else:
            r["file"] = os.path.basename(f)
            reps.append(r)
        if n % 100 == 0:
            print(f"  {n}/{len(man)}", flush=True)
    log(f"parsed {len(reps)} reports, {len(bad)} failed")
    for f, e in bad:
        log(f"  FAILED {f}: {e}")

    best = {}
    for r in reps:
        k = (r["year"], r["week"])
        score = ("corrected" in r["file"].lower(), r["week_end"])
        if k not in best or score > best[k][0]:
            best[k] = (score, r)
    log(f"{len(reps) - len(best)} duplicate (year, week) reports dropped; {len(best)} kept")
    reps = [v[1] for _, v in sorted(best.items())]

    chk = pd.DataFrame([(r["year"], r["week"], int(r["df"]["pools"].sum()), r["t1_pools"]) for r in reps],
                       columns=["year", "week", "t3_pools", "t1_pools"])
    ok = chk.dropna()
    mism = ok[ok.t3_pools != ok.t1_pools]
    log(f"Table 3 pool sum equals Table 1 WNV pools in {len(ok) - len(mism)} of {len(ok)} reports")
    for y, g in mism.groupby("year"):
        log(f"  {int(y)}: {len(g)} mismatched weeks, e.g. " +
            "; ".join(f"wk{int(x.week)} t3={int(x.t3_pools)} t1={int(x.t1_pools)}" for _, x in g.head(6).iterrows()))

    cum = pd.concat([r["df"].assign(year=r["year"], week=r["week"], week_end=r["week_end"])
                     for r in reps if len(r["df"])], ignore_index=True)
    cum = cum.groupby(["year", "week", "week_end", "county"], as_index=False)[FIELDS].sum()
    cum.to_csv(os.path.join(OUT, "tx_wnv_county_cumulative.csv"), index=False)

    weeks = sorted({(r["year"], r["week"]) for r in reps})
    ends = {(r["year"], r["week"]): r["week_end"] for r in reps}
    out = []
    for y in sorted({w[0] for w in weeks}):
        ys = cum[cum.year == y]
        wk_list = sorted(w for (yy, w) in weeks if yy == y)
        counties = sorted(ys.county.unique())
        if not counties:
            continue
        full = (ys.set_index(["week", "county"])[FIELDS]
                .reindex(pd.MultiIndex.from_product([wk_list, counties], names=["week", "county"]), fill_value=0)
                .reset_index().sort_values(["county", "week"]))
        d = full.groupby("county")[FIELDS].diff()
        first = full.groupby("county").cumcount() == 0
        d[first] = full.loc[first, FIELDS]
        neg = int((d < 0).sum().sum())
        if neg:
            log(f"  {y}: {neg} negative weekly differences (report corrections), clipped to 0")
        d = d.clip(lower=0).add_prefix("new_")
        out.append(pd.concat([full[["week", "county"]].assign(year=y), full[FIELDS].add_prefix("cum_"), d], axis=1))
    wk = pd.concat(out, ignore_index=True)
    wk["week_end"] = [ends[(y, w)] for y, w in zip(wk.year, wk.week)]
    wk.to_csv(os.path.join(OUT, "tx_wnv_county_week.csv"), index=False)

    for y in sorted({w[0] for w in weeks}):
        have = {w for (yy, w) in weeks if yy == y}
        missing = [w for w in range(min(have), max(have) + 1) if w not in have]
        if missing:
            log(f"  {y}: missing weeks {missing}")
    log(f"wrote {len(wk)} county-week rows, {wk.county.nunique()} counties")
    with open(os.path.join(OUT, "tx_parse_log.txt"), "w") as f:
        f.write("\n".join(LOG))


if __name__ == "__main__":
    main()
