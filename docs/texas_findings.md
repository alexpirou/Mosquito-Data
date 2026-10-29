# Texas: findings (October 2026)

This covers 2019-2025. Every number comes from `data_tx/tx_results.json`, written by `scripts/tx_analysis.py`. The pipeline is `tx_download_reports.py` -> `tx_parse_reports.py` -> `inat_pull.py --state TX` -> `tx_build.py` -> `tx_analysis.py`.

Texas is the one place where we can test sightings against **WNV-positive mosquito pools**, which turn positive before people get sick and are far more numerous than human cases. DSHS reported 117 new positive pools in 2019 and 3,314 in 2024.

## What the data looks like

- **Sightings:** 11,325 mosquito sightings from 3,545 observers across 254 counties. Only 8.4% are identified as Culex and 45.5% are research grade. 71% of observers posted a single sighting.
- **Sparse:** only 5.3% of county-weeks have any sighting, and 1.6% have a new positive pool.
- **Crowd:** the top 10 users post 16.8% of sightings. The 29 users with 50 or more sightings post 29%. Of the 511 county-weeks with 5 or more sightings, 18% had one person supplying at least 80% of them.
- **Where the pools are:** 46 counties ever reported a positive pool, and four of them (Tarrant, Dallas, Harris, Montgomery) hold most of them. Counties not listed in a report had no activity, so they are zeros, but we can't tell "no positive pools" from "no testing".
- **Seasonality:** sightings peak around MMWR week 18 (center of mass week 29), while pools peak around week 32 (center 33) and human onsets around week 40. Sightings peak in spring, ahead of the pools, but the weekly shapes barely correlate (r = 0.08).

## Pool reports

The DSHS county tables are year-to-date, so weekly new counts are the difference between consecutive reports. We checked the parsed county totals against each report's own summary table, and they match in every week from 2018 on except 2022 week 42 (off by 2 pools). Two limits:

- **Dates are reporting weeks.** A pool counts in the week DSHS reported it, not the week it was collected, so any lead is measured from a date that already lags the mosquitoes.
- **Missing reports.** 2019 week 18 and 2023 week 3 have no report. Outcomes are left missing for any week whose difference would be lumped.

## Main test: do sightings lead pools?

We used 35 counties: those with a positive pool and at least 20,000 iNaturalist observations over 2019-2025. They hold 99.3% of pools and 74% of sightings. The signal is the mosquito share of all iNaturalist observations over a trailing 4-week window, on a log scale. Texas has no insect baseline (as in the national analysis), so the baseline is all observations. The outcome is log(1 + new pools) in MMWR weeks 20-42.

The strictest model gives each county its own usual seasonal shape (county x week-of-year) plus statewide year-week effects, so only a county's deviation from its own normal week counts.

| Signal leads pools by | Effect | p |
|---|---|---|
| 0 weeks | 0.032 | 0.10 |
| 2 weeks | 0.049 | 0.02 |
| 4 weeks | 0.060 | 0.002 |
| 6 weeks | 0.064 | 0.003 |
| 10 weeks | 0.072 | 0.001 |

An unusually mosquito-heavy stretch is followed by more positive pools. One standard deviation of the signal goes with about 0.056 more on the log(1 + pools) scale, which is small: the average county-week has 1.5 new pools, and only 22% of county-weeks in the season have any.

This is different from California and the national analysis, so we checked it hard:

- **Placebo:** pairing each county's signal with a different year's pools gives a null centered near zero. The real estimate is about 2.8 standard deviations above it at 5 weeks (one-sided p = 0.002 in 500 shuffles; at 2 weeks p = 0.024). Shuffling years together keeps the shared regional weather within a year.
- **Not just the county's own pool history:** adding the county's year-to-date pools to the model barely moves it (0.057, p = 0.0003 at 4 weeks). So it isn't only sightings echoing pools the health department already knows about.
- **Not one year:** dropping any single year leaves the estimate between 0.042 and 0.075 (p from 0.001 to 0.048).
- **Same under two other fixed-effect choices** (a county x year version gives 0.035 at 4 weeks, p = 0.003).
- **Both crowd measures agree:** unique observers per observation gives nearly the same curve.

## What the result does not show

- **It is not a particular lead time.** The effect rises to about 0.06 by 3 weeks and then stays flat out to 10 weeks. That looks like "counties and seasons with more mosquito sightings also have more positive pools", not a sharp early warning. We chose 4 weeks for the robustness checks after seeing the curve, and we looked at many lags and variants, so the p-values are optimistic.
- **It is missing where the pools are.** Within the 10 counties with the most pools, the effect is 0.009 (p = 0.83). Dropping those 10 leaves 0.025 (p = 0.07), and dropping just the top 4 leaves 0.040 (p = 0.02). The result comes mostly from mid-sized counties with few pools.
- **Culex alone shows nothing** (0.019, p = 0.35), consistent with only 8.4% of sightings being Culex.
- **Human cases show nothing at any lag** (all effects are 0.011 or smaller and none is significant at 5%; the closest is 10 weeks, p = 0.07). The county tables hide counties with fewer than 5 cases a year, so this test is weak.
- **Season-level checks are not significant.** Using May-June sightings to predict a county's season total, the effect is 0.116 (p = 0.20) over 245 county-years. The statewide year-level rank correlation is 0.61 but with only 7 years (p = 0.15). In a leave-one-year-out forecast, adding sightings to a county's own history cuts the error by only 2.5% (1.042 to 1.016).
- **Fewer counties, tighter sample:** counties with at least 100,000 observations give 0.045 (p = 0.10).

## Bottom line

Texas gives a weak hint that the crowd carries some information about mosquito-testing results, which the national and California analyses did not find. But the evidence is a flat, small association rather than an early warning:

- it is not tied to a lead time,
- it is absent in the highest-activity counties where a forecast would matter,
- it does not show up for human cases,
- and it adds almost nothing to a season forecast.

Positive pools are also not independent of testing effort, and both sightings and trapping may track how many people and agencies are looking. The cleanest reading is that sightings measure human attention to mosquitoes as much as mosquito activity.

## Open items

- A real insect baseline for Texas would be better than all observations (as in the national analysis). It needs another iNaturalist pull.
- Trap-night or testing-effort data would help separate "more mosquitoes" from "more testing".
- The 2026 reports are parsed but not used because the iNaturalist pull stops at 2025.
