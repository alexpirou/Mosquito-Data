# California pilot: findings (October 2026)

This covers 2019-2025. All numbers come from `wnv_inat_analysis.ipynb`.

## What the data looks like

iNaturalist mosquito sightings in California grew from 51 in 2015 to about 1,400 a year in 2023-2025. Before about 2018 there isn't enough data to work with, so we limited the analysis to 2019-2025.

Only 18% of sightings are identified as Culex, the genus our proposal focuses on (mostly *Cx. tarsalis*, *quinquefasciatus* and *pipiens*). Another 15% are the invasive *Aedes aegypti* and *albopictus*, which are common in Los Angeles but aren't the main WNV carriers. Only about 24% of sightings are research grade.

No small group of users dominates the data. There were 2,970 observers across 2015-2025, and the top 5 account for about 12% of sightings. We flagged 8 "heavy posters" with 50 or more sightings in 2019-2025, and together they account for 16%.

Coverage varies a lot by county. Los Angeles has at least one sighting in 83% of weeks, San Diego in 70%, Riverside in 63%, Orange in 43% and Sacramento in 24%. The San Joaquin Valley counties are only at 5-14% each, or 41% when grouped as a region. Unfortunately, the counties with the most WNV are some of the sparsest for sightings.

Sacramento looked promising at first, but its fall 2023 spike in sightings came from one person, who posted 348 of the 351 sightings in that window. It's a mix of *Aedes* and *Culex*, which suggests trap catches rather than casual photos. The raw-count correlation in Sacramento was about 0.30, but it drops to about 0.10-0.14 once we count unique observers or leave out heavy posters. This is the "few highly active users" problem our proposal anticipated.

## Weekly timing

We removed seasonality by comparing each week with the same week in the other years, then checked how well sighting anomalies line up with case anomalies 0 to 8 weeks later.

The best correlations were about 0.10-0.25 for Los Angeles, Orange, Riverside and Sacramento. San Diego and the San Joaquin Valley were close to zero.

To see whether any of this could be chance, we shuffled which year's cases lined up with which year's sightings and compared. Nothing clearly beat chance. The closest was Los Angeles at a 2-week lag (correlation 0.15, p of about 0.07).

## Season-level test

We also asked a simpler question: in years when a county's May-July sightings run above its own normal, does it have more cases in August-November? Across 16 counties and 112 county-years, the correlations were 0.06-0.14 and none were significant. Leaving out invasive *Aedes* did best (0.14, p 0.15). Culex only did worst (-0.03), likely because too few sightings are identified to Culex.

## Spatial clustering

This is the fifth measure from our proposal. We counted a sighting as clustered if a different person reported a mosquito within 2 km and within 7 days either side. We left out obscured observations and those with location accuracy worse than 1 km, which kept 82% of sightings.

Clustering turned out to be rare: 6% of sightings at 1 km, 9% at 2 km and 16% at 5 km. San Diego had at least one cluster in 31 weeks, Los Angeles in 28 and Riverside in 18. Orange, Sacramento and the San Joaquin Valley only had 2-4, which is too few to test.

Where we could test it, there was no clear signal. Los Angeles reached about 0.13 at a 2-week lag (p 0.12 on the shuffle test). San Diego was about 0.20, but it only had 9 cases in the whole period, so that number doesn't mean much. At the season level, years with more clustering didn't have more cases (0.10, p 0.31).

## Where this leaves us

In California from 2019 to 2025, mosquito sightings don't clearly come before human WNV cases. It looks like the main reasons are that sightings are sparse where WNV is worst, invasive *Aedes* make up much of the urban data, a few heavy posters can create false signals, and there aren't enough people reporting close together for clustering to work. This fits the part of our proposal about finding the limits of crowdsourced data for disease surveillance.

## Caveats

CDPH weeks are when a case was reported, not when symptoms started, which pushes cases later than the actual infections.

Weeks 52 and 53 in the CDPH file look like year-end catch-up reporting (205 of 1,478 cases in 2019-2025), so we left them out.

Weekly values tend to follow the week before, so ordinary p-values would make results look stronger than they are. That's why we used the year-shuffle test instead. Picking the best of 9 lags also makes the "best lag" correlations look a bit better than they really are.

Seven seasons is a short record. A weak result here shows the limits of this signal, but it doesn't prove there's no signal at all.
