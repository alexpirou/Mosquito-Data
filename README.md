# Crowdsourced Mosquito Sightings as an Early Signal of West Nile Virus

CS/INFO 5306 project by Alexander Pirouznia, Pratheet Joshi and Sriram Sundararajan.

This project asks whether iNaturalist mosquito sightings can give an early signal of West Nile virus (WNV) activity. Our original proposal is in [`docs/project_proposal.md`](docs/project_proposal.md).

The repo has two analyses that answer the question at different scales:

- **National analysis** ([`national/`](national/)): all 50 states and DC, 2016-2025, using CDC ArboNET neuroinvasive cases by state-year, national month and county. This is the main analysis in the course report. Findings: [`docs/national_findings.md`](docs/national_findings.md).
- **California weekly analysis** (repo root): county-weeks in California, 2019-2025, using the California Department of Public Health's weekly human case counts. This is the finer-grained timing test that the national CDC data can't support. Findings: [`docs/california_pilot_findings.md`](docs/california_pilot_findings.md).

The course report and summary slide are in [`report/`](report/).

## What we found

Both analyses reach the same answer: mosquito sightings don't predict where or when West Nile will hit.

Nationally, there's a loose link across whole seasons, so mosquito-heavy years tend to be worse WNV years. But within a state, an unusually mosquito-heavy spring doesn't predict a worse season (p = 0.25), a crane fly placebo looks no different, and adding sightings to a state's own history makes forecasts slightly worse.

In California, weekly sightings don't clearly come before cases. Correlations are about 0.1 to 0.25 and can't be told apart from chance, and the spatial clustering measure from our proposal is too sparse to test in most counties.

The reasons line up between the two:

- **Culex, the main carrier, is a small share of sightings:** 10% nationally and 18% in California. Most sightings are *Aedes* and *Psorophora*, and Culex sightings peak in the fall, after cases have already peaked.
- **County-weeks are almost always empty:** only 2.6% have any sightings nationally, and the counties with the most WNV in California are some of the sparsest.
- **Raw counts measure iNaturalist's growth and state size,** not mosquitoes. Both analyses adjust for this by using shares of overall activity.
- **The crowd is broad nationally, but one person can still fake a local signal.** Nationally the top ten users post only 5.5% of sightings. In California, though, Sacramento's fall 2023 spike came almost entirely from one user (348 of 351 sightings), and it disappeared once we counted unique observers.

In short, people photograph the mosquitoes that bother them, not the ones that carry the virus.

The two pulls also check each other. For California, the national pull and the California pull share 7,168 of about 7,170 non-casual sightings, and both find 18% are Culex.

## Where things stand

| Part | Status |
|---|---|
| National analysis, 2016-2025 (state-year panel, national months, counties, placebo, forecast test) | Done |
| California weekly analysis, 2019-2025 (lag analysis, season-level test, spatial clustering) | Done |
| Course report and slide | Drafted, in `report/` |
| Texas weekly DSHS reports (mosquito pools and human cases) | Script written, not run yet |
| Texas iNaturalist pull | Script written, not run yet |
| Texas analysis with mosquito pools as the outcome | Not started |

## Repository layout

```
README.md
requirements.txt
docs/
  project_proposal.md
  national_findings.md
  california_pilot_findings.md
report/
  main.tex, references.bib, main.bbl   course report (LaTeX)
  WNV_Project2_Report.pdf              compiled report
  WNV_Project2_Slide.pptx              summary slide
  figures/                             figures used in the report
national/                              national analysis (self-contained; run from this folder)
  README.md
  scripts/pull_data.py, build.py, analysis.py
  data/raw/                            iNaturalist and CDC pulls made 2026-10-26
  data/processed/                      analysis tables
  figures/, results/                   outputs, including results.json with every quoted number
wnv_inat_analysis.ipynb                California analysis (already run, results saved in the notebook)
scripts/
  inat_pull.py                         iNaturalist pull for CA or TX
  tx_download_reports.py               downloads and extracts the Texas DSHS weekly reports
data_ca/                               California data (iNaturalist sightings and baseline, CDPH cases)
data_national/                         CDC national yearly totals, 1999-2025 (background only)
figures/                               California figures
```

## Setup

```
pip install -r requirements.txt
```

### Rerunning the national analysis

Run from the `national/` folder. The raw pulls used for the report are already in `national/data/raw/`, so these two steps reproduce every number without downloading again:

```
cd national
python scripts/build.py
python scripts/analysis.py
```

`python scripts/pull_data.py` downloads everything again (about 45 minutes, resumable). The CDC site returns 403 to Python's requests library, so that step uses curl. See [`national/README.md`](national/README.md) for details.

### Rerunning the California analysis

Run from the repository root. The data is already in `data_ca/`, so you can open `wnv_inat_analysis.ipynb` and run it top to bottom.

To refresh the iNaturalist data:

```
python scripts/inat_pull.py --state CA --start 2015 --end 2025 --wnv data_ca/wnv_human_cases.csv
```

This takes about 15 to 30 minutes. `data_ca/cache/places.json` lets it skip the county lookup.

After the original pull, we reassigned each sighting's county from its coordinates using Census boundaries, because iNaturalist's county matching failed for many counties on the first try. A fresh pull uses iNaturalist's own county places, which agreed with the coordinate method for about 98% of sightings.

### Collecting the Texas data (next step)

These two can run at the same time from the repository root:

```
python scripts/tx_download_reports.py
python scripts/inat_pull.py --state TX --start 2019 --end 2025
```

`tx_download_reports.py` takes about 10 to 20 minutes. It downloads around 450 weekly PDFs into `tx_reports/` (left out of git) and writes `tx_reports_extracted.jsonl` and `tx_reports_manifest.csv`.

`inat_pull.py --state TX` takes about 1 to 1.5 hours because Texas has 254 counties. Output goes to `data_tx/`. If it stops partway, run the same command again and it picks up where it left off.

## Next steps

1. **Update the report with the California results.** The report's limitations section says a true lead test at the state-week level wasn't possible with public data. California's weekly county data makes it possible, so that sentence should change, and a short section on the California timing test, the clustering measure and the Sacramento example would fill the gap. The report also suggests weekly state data as future work, which is what the California and Texas pieces do.
2. **Run the two Texas scripts** and check that `inat_pull.py` prints `version 5 (multi-state)` and finds all 254 counties.
3. **Parse the Texas reports.** We need Table 3, the county-level WNV table (mosquito pools, birds, horses, sentinel chickens, and human fever and neuroinvasive cases). Before trusting the numbers, check them against a couple of the real PDFs:
   - Our understanding is that the county tables are year-to-date totals, so each week's new counts would be the difference between consecutive reports. This still needs to be confirmed.
   - Text extraction mixes up some columns in Table 3, so check the column mapping against the PDF layout.
   - File names and folders change from year to year, some weeks are missing or listed twice, and some reports are "corrected" versions. Use the latest version of each week.
   - Counties with fewer than 5 human cases in a year aren't reported at the county level. Mosquito pool counts shouldn't have this problem.
4. **Run the Texas analysis** with positive mosquito pools as the main outcome. Pools turn positive before people get sick and there are far more of them. It's worth borrowing the national analysis's design here: mosquitoes as a share of insect observations, the crane fly placebo, and the forecast comparison against each county's own history.
5. **Possible extensions:** weather data (PRISM, Daymet or NOAA), which could test whether hot summers drive both the sightings and WNV, and a research-grade Culex-only version.

## Key definitions

**Mosquito share (national):** (mosquito records + 0.5) / (insect records + 1), from the same place and months, analyzed on a log scale.

**Per 1k (California):** per 1,000 total iNaturalist observations in that county-week. The national analysis uses insects as the baseline, which is the better choice for future work.

**Weeks** are CDC (MMWR) weeks, Sunday to Saturday. Week 1 is the first week with at least 4 days in the new year.

**Observers** are the unique users who reported a mosquito in a given place and period, so one very active user can't inflate the number.

**Anomaly (California)** is a week's value minus the average for that same week in the other years.

**Heavy posters (California)** are users with 50 or more mosquito sightings in 2019-2025.

**Clustered sightings (California)** are ones where a different person reported a mosquito within 2 km and within 7 days either side.

## Data sources and notes

**iNaturalist:** the [iNaturalist API](https://api.inaturalist.org/v1/docs/), following its [recommended practices](https://www.inaturalist.org/pages/api+recommended+practices) of about one request per second. The national pull uses verifiable observations only, and the California pull includes all quality grades. Each observation is licensed by the person who posted it, so cite iNaturalist when publishing. These files contain public observation details (numeric user IDs, dates and coordinates) but no photos.

**CDC ArboNET:** [CDC West Nile historic data](https://www.cdc.gov/west-nile-virus/data-maps/historic-data.html), finalized September 16, 2026. County-by-year cases and national cases by onset month, used for the national analysis.

**California human cases:** California Department of Public Health, [West Nile Virus Cases, 2006-present](https://catalog.data.gov/dataset/west-nile-virus-cases-2006-present) (CC-BY). The week is when a case was reported, not when symptoms started, and the file doesn't separate neuroinvasive cases. Only county-weeks with at least one case are listed, so zeros are filled in during the analysis. Weeks 52 and 53 look like year-end catch-up reporting and are left out of the main analysis.

**Texas:** Texas DSHS [weekly arbovirus activity reports](https://www.dshs.texas.gov/idcu/disease/arboviral/westNile/reports/weekly), 2018 onward.

**County population:** 2016 Census estimates via Plotly's datasets (national analysis).

**County boundaries:** U.S. Census boundaries via [us-atlas](https://github.com/topojson/us-atlas).
