# Robustness Classification Method (`classify_rob`)

## Purpose

This note documents the pragmatic rule-based robustness classification used in the reanalysis sensitivity analysis.

## Inputs used by the classifier

For each AMR phenotype, climate variable, and sample definition (`full` or `common`), the classifier combines:

1. Smooth-term support from the GAMM models
   - significance (`Sig`)
   - whether the term was effectively penalized out (`PenOut`)
2. Pairwise curve agreement between the primary source and each alternative source
   - mean curve Pearson correlation (`Curve_r_Mean`)
   - shape agreement (`Shape_Agreement_Rate`)
3. Percentile-based contrast stability
   - direction agreement for the P90 vs P10 contrast (`Direction_Agreement_Rate`)
   - mean ratio of alternative to primary `|Δlog(OR)|` (`Magnitude_Ratio_Mean`)
   - whether a reversal occurs at a non-trivial effect size (`Reversal_Substantive_Rate`)

## Independent alternative sources

- For temperature, precipitation, and wet days: `ERA5` and `MERRA-2` are both treated as alternative sources.
- For relative humidity: only `MERRA-2` is treated as an independent alternative source, because the primary analysis and the ERA5 analysis both use ERA5 humidity.

## Definitions

- `Primary supported`: the primary-source smooth term is significant and not penalized out.
- `Sig_Retained_Rate`: the proportion of independent alternative sources that remain significant and not penalized out.
- `Curve_r_Mean`: the mean Pearson correlation of the smoothed OR curves between Primary and the independent alternative sources.
- `Shape_Agreement_Rate`: the proportion of independent alternative comparisons with the same broad curve class.
- `Direction_Agreement_Rate`: the proportion of independent alternative comparisons with the same sign of the P90–P10 effect contrast.
- `Magnitude_Ratio_Mean`: the mean of `|Δlog(OR)_alt| / |Δlog(OR)_primary|`, calculated only when the primary contrast magnitude exceeds 0.05 in absolute value.
- `Reversal_Substantive_Rate`: the proportion of independent alternative comparisons showing both a direction reversal and a non-trivial percentile-based effect contrast. A reversal is treated as substantive only when the larger of the two `|Δlog(OR)|` values is at least 0.05 and the smaller is at least 0.03.

## Classification rules

### 1. `Robust (consistently NS)`

Assigned when the primary association is not supported and all independent alternative sources are also non-significant or penalized out.

### 2. `Product-sensitive (emergent)`

Assigned when the primary association is not supported but at least one independent alternative source yields a retained significant association.

### 3. `Product-sensitive (lost)`

Assigned when the primary association is supported but all independent alternative sources lose that association.

### 4. `Product-sensitive (reversal)`

Assigned when the primary association is supported and at least one independent alternative source reverses the sign of the percentile-based exposure contrast at a non-trivial effect magnitude.

### 5. `Robust`

Assigned when all of the following are satisfied:

- `Sig_Retained_Rate = 1.00`
- `Curve_r_Mean >= 0.80` when available
- `Shape_Agreement_Rate >= 0.50` when available
- `Magnitude_Ratio_Mean` is between `0.50` and `2.00` when available

### 6. `Robust (minor variation)`

Assigned when all of the following are satisfied:

- `Sig_Retained_Rate >= 0.50`
- `Curve_r_Mean >= 0.60` when available
- `Shape_Agreement_Rate >= 0.33` when available
- `Magnitude_Ratio_Mean` is between `0.25` and `3.00` when available

### 7. `Attenuated`

Assigned to all remaining cases in which the primary association is supported, but the evidence from independent alternative sources is weaker, partially retained, or more variable in magnitude or shape.

## Interpretation guidance

- `Robust` or `Robust (minor variation)` supports retaining the original interpretation.
- `Attenuated` supports keeping the direction of the finding while softening the wording.
- `Product-sensitive (lost)` or `Product-sensitive (reversal)` classifies the finding as source-sensitive.

## Rationale

These thresholds are pragmatic organizing rules rather than formal inferential cutoffs. They are intended to separate clearly stable, moderately stable, attenuated, and source-sensitive findings in a transparent and reproducible way.
