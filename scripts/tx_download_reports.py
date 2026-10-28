"""
Download Texas DSHS weekly arbovirus activity reports (2018 onward) and
extract their text and tables, for the Mosquito Sightings / West Nile project.

What it does
  1. Reads the DSHS weekly reports page and collects every PDF link.
     (File names and folders change from year to year; this handles that by
     taking the links straight from the page.)
  2. Downloads each PDF into tx_reports/<year>/, skipping ones you already have.
  3. Extracts the text and tables from every PDF into one small file,
     tx_reports_extracted.jsonl, which you upload to Claude for parsing.

Usage
  pip install requests pdfplumber
  python tx_download_reports.py

Upload afterwards:
  tx_reports_extracted.jsonl  (a few MB)
  tx_reports_manifest.csv
  and any two PDFs from different years (e.g. one 2019, one 2025) so the
  table layout can be checked against the real thing.

Runtime: roughly 10 to 20 minutes (about 450 PDFs, one per second).
"""

import csv
import json
import os
import re
import time
from urllib.parse import urljoin, quote, urlsplit, urlunsplit

import requests

INDEX = "https://www.dshs.texas.gov/idcu/disease/arboviral/westNile/reports/weekly"
HEADERS = {"User-Agent": "MosquitoSightingsResearch/0.2 (Cornell student project)"}
OUT_DIR = "tx_reports"
PAUSE = 1.0

MONTHS = {m: i for i, m in enumerate(
    ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"], 1)}
NAME_RE = re.compile(r"(\d{1,2})([A-Za-z]+?)(\d{4})[_\- ]?week[_\- ]?(\d{1,2})", re.I)


def clean_url(href):
    url = urljoin(INDEX, href.strip())
    parts = urlsplit(url)
    # Encode spaces and other unsafe characters in the path, keep existing %xx escapes
    return urlunsplit(parts._replace(path=quote(parts.path, safe="/%")))


def parse_name(fname):
    m = NAME_RE.search(fname)
    if not m:
        return None, None, None
    day, mon, year, week = m.groups()
    month = MONTHS.get(mon[:3].lower())
    report_date = f"{year}-{month:02d}-{int(day):02d}" if month else None
    return report_date, int(year), int(week)


def collect_links():
    html = requests.get(INDEX, headers=HEADERS, timeout=60).text
    hrefs = re.findall(r'href="([^"]+?\.pdf)"', html, flags=re.I)
    seen, links = set(), []
    for h in hrefs:
        url = clean_url(h)
        if url.lower() in seen:
            continue
        seen.add(url.lower())
        links.append(url)
    print(f"Found {len(links)} PDF links on the reports page")
    return links


def download(links):
    rows = []
    for i, url in enumerate(links, 1):
        fname = os.path.basename(urlsplit(url).path).replace("%20", " ")
        report_date, year, week = parse_name(fname)
        folder = os.path.join(OUT_DIR, str(year or "unknown"))
        os.makedirs(folder, exist_ok=True)
        path = os.path.join(folder, fname)
        status = "exists"
        if not os.path.exists(path):
            time.sleep(PAUSE)
            try:
                r = requests.get(url, headers=HEADERS, timeout=120)
                if r.status_code == 200 and r.content[:4] == b"%PDF":
                    with open(path, "wb") as f:
                        f.write(r.content)
                    status = "downloaded"
                else:
                    status = f"failed HTTP {r.status_code}"
            except requests.RequestException as e:
                status = f"failed {type(e).__name__}"
        rows.append({"url": url, "file": path, "report_date": report_date,
                     "year_in_name": year, "week_in_name": week, "status": status})
        if i % 25 == 0 or status.startswith("failed"):
            print(f"  [{i}/{len(links)}] {fname}: {status}")
    with open("tx_reports_manifest.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    ok = sum(r["status"] in ("downloaded", "exists") for r in rows)
    print(f"Downloaded or already had {ok} of {len(rows)} PDFs. See tx_reports_manifest.csv")
    return rows


def extract(rows):
    import pdfplumber
    n = 0
    with open("tx_reports_extracted.jsonl", "w", encoding="utf-8") as out:
        for r in rows:
            if not os.path.exists(r["file"]):
                continue
            try:
                with pdfplumber.open(r["file"]) as pdf:
                    pages = []
                    for pno, page in enumerate(pdf.pages, 1):
                        pages.append({
                            "page": pno,
                            "text": page.extract_text() or "",
                            "tables": page.extract_tables(),
                        })
            except Exception as e:
                pages = [{"page": 0, "error": f"{type(e).__name__}: {e}"}]
            out.write(json.dumps({"file": os.path.basename(r["file"]),
                                  "report_date": r["report_date"],
                                  "year_in_name": r["year_in_name"],
                                  "week_in_name": r["week_in_name"],
                                  "pages": pages}) + "\n")
            n += 1
            if n % 50 == 0:
                print(f"  extracted {n}")
    size = os.path.getsize("tx_reports_extracted.jsonl") / 1e6
    print(f"Extracted {n} reports into tx_reports_extracted.jsonl ({size:.1f} MB)")


if __name__ == "__main__":
    links = collect_links()
    rows = download(links)
    extract(rows)
