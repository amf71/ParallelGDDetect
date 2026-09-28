# ParallelGDDetect
A package to detect distinct genome doubling events which have occurred in different samples from the same tumour.

A function is provided which takes as input the number of genome doublings (GDs) for each region as estimated from the genome
wide copy number and mutation copy numbers in each region and devolves which GDs across samples are part of
the same event and which present distinct events as indicated by the doubling of subclonal mutations which
will occur when the subclonal mutation arises before a given subclonal GD event. We recommend plotting the mutation
copy numbers alongside the allele specific copy number across the genome to verify these calls as well as carefully
checking the ploidy solutions determined for each sample. Default thresholds have been set for best performance
on TRACER NSCLC exome sequencing data. Particularly for whole genome sequencing data, thresholds
may need modification.

## Installation & loading

You can use devtools::install_github() to install cloneMap from this repository:
Temporary authenication token for reviewers included

```R
devtools::install_github('amf71/ParallelGDDetect')
```

load package:

```R
library(ParallelGDDetect)
```

## Inputs

input A input table describing each mutation in a tumour or set of tumours with the following columns:
* tumour_id: A unique identifier for each tumour
* chromosome: The chromosome in which the mutation is present
* position: the base position of the mutations in the chromosome
* ref: The reference base
* alt: the alternate/variant base
* cluster_id: A unique id for each mutational cluster representing past or current subclones of the tumour. This can be calculated using tools like PyClone. While this has not been tested, it may also be effective for this method to simply use a clustering of mutations based on presence or absence in each region while limiting the input to mutations with at least 0.75 mut_cpn in one region.
* is_clonal_cluster: Is the cluster the clonal cluster
* sample_id: A unique identifier for each sample in each tumour
* mut_cpn: The estimated non-integer mutation copy number for each mutation. This is calculated by tools like PyClone where joint inferences of CCF and multiplicity are made during mutation clustering or can be calculated more simply using the equation on page 14 of supplementary appendix 1 of Jamal-Hanjani et al 2017 NEJM. AS this methods repliesonly on mutations which are clonal in a given region (even if they are subclona accross the tumour as a whole)
either method will be appropriate.
* MajCN: The copy number of the major allele at the mutated locus.
* num_gds: The number of genome doublings that have occured in each sample as estimated from the ploidy

We recommend using using thresholds of >= 50% of the genome with at least Major allele copy number >= 2 for determination of whether a First GD (from ploidy 2 to 4) has occured and a thshold of >= 50% of the genome with at least Major allele copy number >= 3 for determination of whether a seocnd GD (from ploidy 3-4 to 6-8) has occured. These thesholds account for the higher frequency of losses rather than gains the are known to occur after a GD event and are similar to those used in carter et al 2012 Nature biotechnology which first published the ABSOLUTE tool. 

# Example Run
## Run on example data (loaded with package)
output <- detect_par_gd( example_data )

# Outputs

A list is returned with three objects which each describe the genome doubling events over all the tumours that were inputted:
* GDs_per_tumour: Description of genome doubling events for each tumour (one row per tumour)
  * First_GD: Is there a first GD event anywhere in the tumour (ie an event that would modify ploidy from 2 to 4)
  * Second_GD: Is there a second GD event anywhere in the tumour (ie an event that would modify ploidy from 3-4 to 6-8)
  * num_first_gd: number of first genome doublings across the tumour (clonal or subclonal)
  * num_second_gd: number of second genome doublings across the tumour (clonal or subclonal)
  * num_clonal_gd: number of clonal genome doublings across the tumour (first or second)
  * num_subclonal_gd: number of subclonal genome doublings across the tumour (first or second)
  * First_GD_homogen: Whether all regions have any first GD event (from ploidy 2 to 4) 
  * Second_GD_homogen: Whether all regions have any second GD event (from ploidy 3-4 to 6-8) 
  * GD_status_homogen: Whether all regions have had the same number of GD events (ie will have ~ the same ploidy)
  * GD_statuses: A common separated list of all the GD states present in the tumour (from 0, 1 or 2)
  * frac_0_gd_regions: The fraction of regions with no GD event
  * frac_1_gd_regions:  The fraction of regions with a first GD event (from ploidy 2 to 4)
  * frac_2_gd_regions:   The fraction of regions with a second GD event (from ploidy 3-4 to 6-8)
   
* GDs_per_region: Description of genome doubling events for each region (one row per region)
  * gd_clusters: A common separated list of the subclonal clusters in a given region which were identified with enough mutations at copy number 2 that some of the mutations probably occurred before a subclonal GD event
  * gd_events: A common separated list of the GD event IDs present in a given region
  * First_GD: The GD event ID for the first GD (from ploidy 2 to 4) in a given region or 'No GD' if there was no first GD 
  * Second_GD: The GD event ID for the second GD (from ploidy 3-4 to 6-8) in a given region or 'No GD' if there was no second GD 

* GDs_events: Description of genome doubling events for each tumour (one rwo per event)
  * GD_event: The event ID unique within each tumour
  * GD_event_id: The event ID concatonated with the tumour id hence is unique accross a cohort
  * tumour_id: A id for the tumour
  * is_clonal: Whether the event is clonal (present in all regions) or subclonal (pesent in a subset of regions)
  * num_regions_gd_event: Number of regions in which the GD is present
  * total_regions: Total number of regions in the tumour which the GD event is present
  * is_subclonal_mutation_supported: Whether there is a doubled subclonal mutation cluster supporting a subclonal GD (TRUE) or if this subclonal GD is inferred only from the ploidy (FALSE)
  * clusters: If there are supporting subclonal mutation clusters which are they (common seperated list, otherwise NA if none)

## Bootstrap stability filtering

By default, `detect_par_gd()` uses the discovery and check thresholds directly, without bootstrapping. Set `bootstrap_discovery = TRUE` to assess how consistently each discovery call is reproduced when mutations are resampled with replacement within each tumour, sample, and cluster. Only tumour/sample/cluster combinations that pass the discovery thresholds on the original data are evaluated for bootstrap stability.

```R
output <- detect_par_gd(
  example_data,
  bootstrap_discovery = TRUE,
  n_boot = 500,
  seed = 1
)
```

Bootstrap filtering builds the standard result tables in two passes:

1. **Stable anchors:** A tumour/sample/cluster must pass the discovery thresholds and the selected stability test to be called from its own evidence.
2. **Check-threshold second chance:** If a cluster has at least one stable anchor anywhere in the tumour, other samples of that cluster can be called using the lower check thresholds, without needing to pass the stability test themselves. Anchors retain their validated calls and are not re-evaluated against the check thresholds. Without a stable anchor, no sample of the cluster can be called through this second pass.

Thus, a stable call in one region can support calls in other regions of the same cluster, even when those regions do not independently pass discovery. Returned `GDs_per_tumour`, `GDs_per_region`, `GDs_events`, `mut_counts`, and `mut_counts_not_called` are built from this filtered call set.

The bootstrap and stability arguments are:

| Argument | Default | Description |
| --- | --- | --- |
| `bootstrap_discovery` | `FALSE` | Turn bootstrap stability filtering on or off. |
| `n_boot` | `500` | Number of resamples per tumour/sample/cluster. |
| `seed` | `1` | Random seed for reproducibility. |
| `stability_threshold` | `0.7` | Minimum stability tested; stability is the fraction of resamples reproducing the discovery call. |
| `stability_alpha` | `0.05` | Significance level for the stability rule. With the binomial test enabled, this is the adjusted q-value cutoff. |
| `stability_use_binom_test` | `TRUE` | Use a one-sided binomial test against the null that true stability is at most `stability_threshold`, with p-values adjusted by the selected method. A call is stable when `q_value < stability_alpha`. If `FALSE`, use the Wilson interval rule: its lower bound must exceed `stability_threshold`. |
| `stability_adjust_method` | `"BH"` | P-value adjustment method passed to `stats::p.adjust` when using the binomial test. |

When bootstrapping is enabled, the returned list also includes:

* `bootstrap_discovery`: One row per original-data discovery call, with `k` (resamples reproducing the call), `B` (resamples run), `stability` (`k / B`), Wilson interval bounds (`ci_low`, `ci_high`), and the binomial test fields (`p_value`, `q_value`) when applicable. `keep_call` records whether `ci_low > stability_threshold`; `is_stable` records the rule actually selected by `stability_use_binom_test`.
* `bootstrap_discovery_settings`: The parameters and settings used for bootstrapping, including the resampling unit and the requirement for a stable cluster anchor before applying the check thresholds.
* `unstable_calls_GDs_per_tumour`, `unstable_calls_GDs_per_region`, `unstable_calls_GDs_events`, `unstable_calls_mut_counts`, and `unstable_calls_mut_counts_not_called`: The calls from the unfiltered discovery-plus-check pipeline that were excluded from the filtered results. These include calls lost to stability filtering and check-phase calls that had no stable anchor. They can be used to inspect calls rejected by the filter.

## Reproducing original publication

In silico benchmarking presented in (Frankell et al., 2023)[https://doi.org/10.1038/s41586-023-05783-5] has been performed with the previous version of this programme. In order to reproduce it, please use 1af16a8355d594cc6cb91cfa3042c6292f5631ae commit has, i.e.'
```bash 
git checkout -b original-publication 1af16a8355d594cc6cb91cfa3042c6292f5631ae
```
