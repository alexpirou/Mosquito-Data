# National analysis: findings (October 2026)

This covers all 50 states and DC, 2016-2025. The code, data and every quoted number are in [`national/`](../national/), and the full write-up is the course report in [`report/`](../report/).

## Data

We pulled every verifiable US mosquito (Culicidae) observation from iNaturalist for 2015-2025, which came to 76,129 observations from 23,144 people. For a baseline, we also pulled monthly counts of all insect observations per state (about 30 million in total), crane fly observations, and yearly insect counts for the 2,357 counties with at least one mosquito record.

For West Nile virus we used the CDC's ArboNET historic data (finalized September 16, 2026), which has cases by county and year and national cases by month of onset. Our outcome is neuroinvasive cases, since the CDC says milder cases are under-reported unevenly and shouldn't be compared across places or years.

Sightings were placed in counties by coordinates against Census boundaries (99.5% placed). For the 11.8% with hidden locations, we used iNaturalist's own county. State totals came within 0.3% of iNaturalist's counts, and county case totals added up exactly to the CDC's national totals.

## Why the design changed from the proposal

Only 10.0% of mosquito observations are Culex. Most are *Aedes* (39.6%) or *Psorophora* (11.9%), and 11.7% are never identified past the family level. So all mosquitoes became the main signal, with Culex as a sub-signal.

County-weeks are almost always empty. The median county has only 6 mosquito observations in ten years, and only 2.6% of county-weeks have any. The CDC's public data also only goes down to county-year. We moved to state-years, national months and a pooled county comparison.

Raw counts mostly measure the app. Mosquito observations grew about fourteen-fold from 2016 to 2025 while WNV showed no trend. We divided each mosquito count by insect observations from the same place and months (the "mosquito share"), and compared each state with its own normal level and with the national year.

## Main question and method

Within a state, does an April-July mosquito share that's higher than usual, in a year that isn't unusual nationally, go along with more neuroinvasive cases that year? April-July comes before the August-September peak in cases.

We fit a regression of log(1 + neuroinvasive cases) on each April-July signal, with state and year fixed effects and standard errors clustered by state (508 state-years). Crane flies, a mosquito look-alike that carries no disease, serve as a placebo.

## Results

**Seasonal timing.** Mosquito observations center on mid-July, about 1.1 months before neuroinvasive onsets center in mid-August. Most of that lead is just when people photograph insects in general. Culex observations peak in October, after the August peak in cases.

**Year to year nationally, there's a loose link.** The April-July mosquito share and that year's cases had a rank correlation of 0.35 (p = 0.33, n = 10), and July-September did better at 0.58 (p = 0.08). With monthly anomalies, an above-normal mosquito share was followed by more cases one and two months later (r = 0.39 and 0.44, permutation p = 0.02). Once each year's average was also removed, that dropped to 0.05-0.15 and wasn't significant, so the link is about whole seasons being mosquito-heavy, not month-to-month timing.

**Within states, there's no signal.** Every mosquito signal had a positive coefficient and the crane fly placebo a negative one, but none were significant and none explained even 1% of the remaining variation. Counting observers or narrowing to Culex didn't help.

| April-July signal | Coefficient | SE | p | Within-R2 |
|---|---|---|---|---|
| Mosquito share of insect observations | 0.072 | 0.063 | 0.25 | 0.003 |
| Mosquito observers per insect observation | 0.058 | 0.089 | 0.51 | 0.001 |
| Culex share of insect observations | 0.050 | 0.038 | 0.19 | 0.003 |
| Raw mosquito observations (log) | 0.090 | 0.073 | 0.22 | 0.005 |
| Unique mosquito observers (log) | 0.098 | 0.106 | 0.35 | 0.003 |
| Placebo: crane fly share | -0.058 | 0.086 | 0.50 | 0.001 |
| Mosquito share, full year (concurrent) | -0.005 | 0.065 | 0.93 | 0.000 |

**Robustness checks told the same story.** Restricting to 2020-2025 gave the largest estimate (0.25, p = 0.08). A negative binomial model looked significant until standard errors were clustered by state (p = 0.10). Dropping year effects made the signal significant (0.13, p = 0.03), which matches the national pattern. 28 of 46 states had a positive within-state correlation, not much better than a coin flip (p = 0.18).

**As a forecast, there's no gain.** Adding the mosquito anomaly to each state's own average made leave-one-year-out forecasts slightly worse (mean absolute error 0.687 vs. 0.682).

**The raw-count trap.** Pooled across state-years, raw mosquito counts correlate with cases at 0.23 (p < 0.001), but that's state size. As a share, it drops to -0.04.

**Counties.** Counties with a higher mosquito share have slightly higher average incidence (0.066, p = 0.002), but within their own state, controlling for population and insect activity, the effect disappears (0.006, p = 0.57).

**The crowd itself isn't the problem nationally.** The 75,798 observations from 2016-2025 came from 23,078 people. 69.4% posted only one mosquito, and the top ten users account for just 5.5%.

## Where this leaves us

The signal is real but weak, and it only shows up nationally, across whole seasons. People mostly photograph *Aedes* and *Psorophora*, and the Culex that carry the virus are a tenth of sightings and get posted most in the fall. A hot summer could plausibly produce both more noticeable mosquitoes and more WNV at the same time, which would explain why the link appears across seasons but not within states. Testing that would need weather data.

## Limitations

Ten years means only ten national seasons. State shares before about 2019 are mostly noise. Public CDC data only has onset timing nationally and yearly totals by state and county. Cases are counted by county of residence, not where the person was bitten. 60% of observations aren't research grade. County boundaries were generalized, and county populations come from a single year.
