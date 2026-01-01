[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18116005.svg)](https://doi.org/10.5281/zenodo.18116005)

# Two-Sample Permutation Tests with Arbitrary Statistics for Paired and Unpaired Designs in MATLAB

Two-sample permutation (randomisation) tests for MATLAB, implementing label-flip permutation for independent samples and sign-flip permutation for paired or repeated measures designs, with support for arbitrary user-defined test statistics.

---

## What this is

This repository provides two small, focused MATLAB functions for **two-sample inference**:
- `labelFlip_permTest2_unpaired.m` — two-sample permutation test for **independent samples**
- `signFlip_permTest2_paired.m` — two-sample permutation test for **paired / repeated-measures samples**

The functions are intentionally with no plotting, UI or fixed statistics. All examples and usage details are documented directly in the function headers.

---

## What problem this solves

Most existing MATLAB permutation-test utilities are:
- hard-coded to **difference in means**
- limited to unpaired designs
- difficult to extend to non-standard scientific questions

This implementation treats permutation testing as a **general two-sample inference engine**. You supply:
- two samples (already aggregated to the correct experimental unit)
- a statistic that directly answers *your scientific question*

The permutation test supplies:
- a valid null distribution
- a p-value with exact Type I error control under exchangeability

---

## Why arbitrary statistics matter

Permutation tests do not require the statistic to be a mean. Any scalar-valued function comparing two samples that respects exchangeability under the null can be tested, including:
- **Measures of variability or heterogeneity** (e.g. variance, standard deviation, MAD, interquartile range)
- **Threshold-based summaries** (e.g. proportion of observations exceeding a predefined cutoff)
- **Nonlinear distributional summaries** (e.g. medians, quantiles, ratios)
- **Robust location or scale estimators** (e.g. trimmed means, winsorized means)
- **Domain-specific scalar metrics** derived from the sample distributions

This enables questions such as:
> “Does condition A increase variability or heterogeneity, even if the mean is unchanged?”

to be tested directly, using the same inference machinery.

---

## Paired vs unpaired designs

The two two-sample designs are implemented **explicitly and separately**:
- **Unpaired**: label permutation across independent samples  
- **Paired**: sign-flip permutation of within-pair differences

This avoids a common error where paired data are analyzed using unpaired logic, violating exchangeability under the null.
For paired sign-flip tests: statistics that change sign when within-pair differences are reversed (e.g. mean or median) can be tested using either one- or two-sided tests. Statistics that are unchanged by sign reversal (e.g. variance or MAD) are typically tested one-sided.

---

## Statistical assumptions

The only assumption is:
> **Exchangeability of the input observations under the null hypothesis.**

This requires that each input value corresponds to an **independent experimental unit** (e.g. per-animal summary statistics).

Hierarchical or nested designs are intentionally not handled automatically. The users must preprocess their data to proxide values at the correct experimental unit before applying the permutation test.

---

## References

- Fisher, R. A. (1971). The Design of Experiments (9th ed.). New York: Hafner Press. 
- Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be zero: Calculating exact P-values when permutations are randomly drawn. *Stat. Appl. Genet. Mol. Biol.*, 9(1), Article 39. https://doi.org/10.2202/1544-6115.1585
- Pitman, E. J. G. (1937). Significance tests which may be applied to samples from any populations. *Supplement to the Journal of the Royal Statistical Society, Series B (Statistical Methodology)*, 4(1), 119–130.

---

## License

MIT License.

---

## Author

Gokul Rajan.

---

## Citation

Gokul Rajan. (2026). labelFlip and signFlip: Two-Sample Permutation Tests in MATLAB with Arbitrary Statistics for Paired and Unpaired Designs (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.18116005
