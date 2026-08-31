
Bootstrap Discovery Stability in `detect_par_gd`

This note documents how the optional bootstrap in `detect_par_gd()` works and how to interpret its outputs.

## What this bootstrap does

When `bootstrap_discovery = TRUE`, the function evaluates the *stability* of mutation-based **discovery-phase** subclonal GD calls (`is_subcl_gd_discovery`).

Important scope:

- It bootstraps only the **discovery phase** thresholds.
- The default (filtered) tables are built through a two-pass procedure instead of the original single check-phase pass:
  1. **Pass 1 (anchors)**: a tumour/sample/cluster must pass the discovery thresholds AND be stable under
     bootstrapping to be called on its own evidence.
  2. **Pass 2 (check, second chance)**: once a cluster has at least one stable anchor somewhere in the tumour,
     every *other* sample of that same cluster gets a second chance at the lower "check" thresholds - without
     needing to pass the stability test itself. If no sample of the cluster is a stable anchor anywhere in the
     tumour, none of its other samples get checked either. Anchors keep their own already-validated call (they
     are not re-evaluated against the check thresholds).
- It splits results into two sets: default tables, where only the filtered (anchor + check-rescued) calls are
  included, and 'unstable_calls_...' tables, which contain whatever the un-filtered discover+check pipeline would
  have called minus what's in the default tables - ie calls (including check-phase rescues with no stable anchor
  backing the cluster) that are dropped specifically because of stability filtering.

## Inputs used by the bootstrap

The bootstrap uses the same original mutation-level input columns required by discovery logic, especially:

- `tumour_id`, `sample_id`, `cluster_id`
- `mut_cpn`, `MajCN`, `num_gds`

It also uses these function arguments:

- `bootstrap_discovery` (enable/disable)
- `n_boot` (number of bootstrap replicates)
- `seed` (reproducibility)
- Discovery thresholds:
  - `discover_mut_cpn_2_threshold`
  - `discover_num_muts_threshold`
  - `discover_frac_2_cpn_muts_threshold`
  - `discover_num_2_cpn_muts_threshold`
- Stability inference settings:
  - `stability_threshold`
  - `stability_alpha`
  - `stability_use_binom_test`
  - `stability_adjust_method`

## Defaults
Bootstrapping is controled via boolean 'bootstrap_discovery' parameter.
n_boot is set to 500 by default.
Stability threshold is set to 0.7, and alpha is set to 0.05 by default. 

## Bootstrap algorithm

For each baseline call (a tumour/sample/cluster that was discovery-positive originally):

1. Run `n_boot` resamples.
2. In each resample, recompute discovery call status.
3. Count how often it remains discovery-positive:

$$
k = \sum_{b=1}^{B} I\{\text{is_subcl_gd_discovery}^{(b)} = TRUE\}
$$

4. Compute empirical stability:

$$
\hat{p} = \frac{k}{B}
$$

5. Compute Wilson confidence interval at level `1 - stability_alpha`.
6. Mark `keep_call = (ci_low > stability_threshold)`.

Optional (`stability_use_binom_test = TRUE`):

- Compute one-sided binomial p-value against `stability_threshold`:

$$
p\_value = P(X \ge k \mid X \sim Binomial(B, p_0 = stability\_threshold))
$$

- Adjust to `q_value` using `stability_adjust_method` (default `BH`).

If  `stability_use_binom_test = TRUE`, results are filtered using q_value and alpha, otherwise they are filtered using `keep_call = (ci_low > stability_threshold)`

## Output tables and fields

When bootstrapping is enabled, the default `GDs_per_tumour`/`GDs_per_region`/`GDs_events`/`mut_counts`/
`mut_counts_not_called` tables carry only the filtered (anchor + check-rescued) calls described above, and several
extra list entries are returned:

1. `unstable_calls_GDs_per_tumour`, `unstable_calls_GDs_per_region`, `unstable_calls_GDs_events`,
   `unstable_calls_mut_counts`, `unstable_calls_mut_counts_not_called`

- Same shape as the default tables, but built from whatever the un-filtered discover+check pipeline would have
  called, minus what's in the default tables.

2. `bootstrap_discovery`

- One row per original discovery-positive tumour/sample/cluster.
- Columns:
  - `tumour_id`, `sample_id`, `cluster_id`
  - Original discovery metrics: `num_muts`, `perc_cn2`, `num_cn2`
  - `k`: number of bootstrap replicates where discovery call is positive
  - `B`: number of replicates (`n_boot`)
  - `stability`: `k / B`
  - `ci_low`, `ci_high`: Wilson CI bounds for stability
  - `keep_call`: TRUE if `ci_low > stability_threshold`
  - Optional: `p_value`, `q_value` (if `stability_use_binom_test = TRUE`)
  - `is_stable`: the rule actually used to determine Pass-1 anchors - `q_value < stability_alpha` if
    `stability_use_binom_test = TRUE`, otherwise `keep_call`

3. `bootstrap_discovery_settings`

- One-row table recording:
  - `n_boot`, `seed`
  - `stability_threshold`, `stability_alpha`, `stability_use_binom_test`, `stability_adjust_method`
  - `bootstrap_unit` (currently `tumour_id x sample_id x cluster_id`)
  - `phase` (currently `discovery_only`)
  - `check_phase_requires_stable_anchor_in_cluster` (currently always `TRUE`)

## How to interpret results

Recommended interpretation for each row:

- High confidence stable call:
  - `stability` high (close to 1)
  - `ci_low` above threshold
  - `keep_call = TRUE`
- Borderline call:
  - `stability` near threshold
  - CI overlaps threshold
  - `keep_call = FALSE`
- Unstable call:
  - low `stability`
  - low `ci_low`
  - usually `keep_call = FALSE`

Practical reading of key fields:

- `stability` is the empirical reproducibility under resampling.
- `ci_low` is the conservative lower-bound estimate of that reproducibility.
- `keep_call` is a strict rule: the lower bound must exceed the threshold.
- Small `p_value`/`q_value` supports stability greater than the threshold under the binomial model.
