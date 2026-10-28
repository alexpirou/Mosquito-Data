# iNaturalist mosquito sightings vs. West Nile virus (CS/INFO 5306 Project 2)

Run from this folder with Python 3.11+ (`pip install pandas numpy scipy statsmodels shapely matplotlib requests`):

1. `python scripts/pull_data.py`   downloads everything into `data/raw/` (about 45 min, mostly iNaturalist's rate limit; resumable).
   cdc.gov serves its CSVs to curl but returns 403 to python-requests, so the CDC step shells out to curl.
2. `python scripts/build.py`       builds the tables in `data/processed/`
3. `python scripts/analysis.py`    writes `figures/*.png`, `results/results.json` and `results/tables.md` (about 2 min; the permutation tests are the slow part)

The raw pulls used for the report (made 2026-10-26) are already in `data/raw/`, so steps 2-3 reproduce the numbers without re-downloading.

## Files

- `data/raw/inat/culicidae_us.jsonl`: 76,129 verifiable US mosquito (Culicidae) observations, 2015-2025
- `data/raw/inat/state_month_hist.json`: monthly counts per state for insects, mosquitoes, Culex, crane flies, all taxa
- `data/raw/inat/county_year_insecta.json`: yearly insect counts per iNaturalist county place
- `data/raw/cdc/*.csv`: CDC ArboNET historic WNV files (finalized 2026-09-16)
- `data/processed/state_year.csv`: the main analysis table (one row per state-year, with Apr-Jul and full-year crowd signals and CDC cases)
- `data/processed/national_month.csv`, `state_month.csv`, `county.csv`, `mosquito_obs.csv`: the other analysis tables
- `results/results.json`: every number quoted in the write-up

## Definitions

- **Mosquito share** = (mosquito records + 0.5) / (insect records + 1), from the same place and months. Analyses use its log.
- **Early season** = April-July, before the August-September peak in human onsets.
- **Outcome** = neuroinvasive WNV cases (CDC advises against comparing milder cases across places or years).
- **Main model**: OLS of log(1 + neuroinvasive cases) on the signal, with state and year fixed effects and state-clustered SEs. Only state-years with at least 200 April-July insect records are kept (508 of 510).
