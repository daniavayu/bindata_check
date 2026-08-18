# Malawi 1997: preliminary top-tail results

## Main finding

The evidence points to two exceptionally influential upper-tail observations,
not to a general problem affecting all high-welfare observations.

- The LIS rule on log welfare with `k = 3` flags 59 observations, representing
  0.435% of the weighted population.
- A Pareto distribution fitted to the top 10% broadly describes most of that
  upper tail.
- The two largest values are far more extreme than the fitted Pareto tail would
  predict. The maximum is about 713 times as frequent as the p90 Pareto model
  predicts; the second-largest value is about 96 times as frequent.

## Influence on reported measures

Removing only those two values reduces the weighted mean from 5.78 to 3.80
2021 PPP USD per person per day (-34.3%). It reduces the Gini from 0.658 to
0.479 (-27.1%), MLD from 0.811 to 0.392 (-51.7%), and Atkinson(2) from 0.647
to 0.463 (-28.5%).

The effect on headcount, FGT(1), FGT(2), and Watts poverty measures is
negligible because the two observations are in the upper tail, not among the
poor.

## Dagum sensitivity check

Dagum Type I was fitted to the entire positive distribution with survey
weights. The global weighted KS distance changes only from 0.01650 with all
records to 0.01641 without the two largest observations. Therefore, the Dagum
fit to the whole distribution is largely governed by the bulk of the data and
does not itself provide evidence to remove the two values.

## Interpretation

These diagnostics do not prove data error. The two records should be checked
against the source survey and welfare-construction process. If they cannot be
validated, this evidence supports a transparent top-tail treatment as a
sensitivity analysis. A decision on top-coding should report both the rule and
its effect on inequality measures.
