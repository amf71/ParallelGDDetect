#====================================================#
#     Function to produce test for Parrallel GDs     #
#    across different samples from the same tumour   #
#====================================================#



# input = example_data
# input = eg_input_all
# 
# discover_mut_cpn_2_threshold = 1.5; check_mut_cpn_2_threshold = 1.25; discover_num_muts_threshold = 20;
# discover_num_2_cpn_muts_threshold = 10; discover_frac_2_cpn_muts_threshold = 0.25;
# check_num_muts_threshold = 10; check_num_2_cpn_muts_threshold = 5; 
# check_frac_2_cpn_muts_threshold = 0.1

#' Detect Parallel subclonal genome doubling
#'
#' Function which takes as input the number of genome doublings (GDs) for each region as estimated from the genome
#' wide copy number and mutation copy numbers in each region and devolves which GDs across samples are part of
#' the same event and which present distinct events as indicated by the doubling of subclonal mutations which
#' will occur when the subclonal mutation arises before a given subclonal GD event. We recommend plotting the mutation
#' copy numbers alongside the allele specific copy number across the genome to verify these calls as well as carefully
#' checking the ploidy solutions determined for each sample. Default thresholds have been set for best performance
#' on TRACER NSCLC exome sequencing data. Particularly for whole genome sequencing data, thresholds, particularly for mutation
#' numbers may need modification.
#'
#' @param input A input table describing each mutation in a tumour or set of tumours with the following columns:
#'     * tumour_id: A unique identifier for each tumour
#'     * chromosome: The chromosome in which the mutation is present
#'     * position: the base position of the mutations in the chromosome
#'     * ref: The reference base
#'     * alt: the alternate/variant base
#'     * cluster_id: A unique id for each mutational cluster representing past or current subclones of the tumour
#'       This can be calculated using tools like PyClone. While this has not been tested, it may also be effective
#'       for this method to simply use a clustering of mutations based on presence or absence in each region while 
#'       limiting the input to mutations with at least 0.75 mut_cpn in one region.
#'     * is_clonal_cluster: Is the cluster the clonal cluster
#'     * sample_id: A unique identifier for each sample in each tumour
#'     * mut_cpn: The estimated non-integer mutation copy number for each mutation. This is calculated by tools like PyClone where
#'       joint inferences of CCF and multiplicity are made during mutation clustering or can be calculated more simply 
#'       using the equation on page 14 of supplementary appendix 1 of Jamal-Hanjani et al 2017 NEJM. AS this methods replies
#'       only on mutations which are clonal in a given region (even if they are subclona accross the tumour as a whole)
#'       either method will be appropriate.
#'     * MajCN: The copy number of the major allele at the mutated locus.
#'     * num_gds: The number of genome doublings that have occured in each sample as estimated from the ploidy
#'       We recommend using using thresholds of >= 50% of the genome with at least Major allele copy number
#'       >= 2 for determination of whether a First GD (from ploidy 2 to 4) has occured and a thshold of >= 50% 
#'       of the genome with at least Major allele copy number >= 3 for determination of whether a seocnd GD 
#'       (from ploidy 3-4 to 6-8) has occured. These thesholds account for the higher frequency of losses
#'       rather than gains the are known to occur after a GD event and are similar to those used in carter et al
#'       2012 Nature biotechnology which first published the ABSOLUTE tool. 
#' 
#' @param discover_mut_cpn_2_threshold The mutation copy number (mut_cpn) threshold to consider a mutation most likely at 2 copies,
#' used when initially discovering whether a cluster shows evidence of a subclonal GD ('discover_' criteria). Default is 1.5.
#'
#' @param check_mut_cpn_2_threshold A lower mutation copy number (mut_cpn) threshold to consider a mutation most likely at 2 copies,
#' used for the 'check_' criteria applied to a cluster in regions where it has already been identified as doubled in at least one
#' other region. If not supplied, defaults to the value of discover_mut_cpn_2_threshold.
#'
#' @param discover_num_muts_threshold The number of mutations required to assess a clone for evidence of a subclonal GD events using 
#' the mutation copy number (mut_cpn). Default is 20. 
#' 
#' @param discover_frac_2_cpn_muts_threshold The fraction of 2 mutation copy number mutations required to consider that some mutations in a 
#' cluster have most likely occurred before a subclonal genome doubling event in a region. A lower threshold is later applied for other regions
#' where a cluster has been called doubled in at least 1 region.  Default is 0.25 
#' 
#' @param check_frac_2_cpn_muts_threshold A lower threshold for the fraction of 2 mutation copy number mutations required to consider that some mutations in a 
#' cluster have most likely occurred before a subclonal genome doubling event in a region if the cluster has already been identified as doubled in other regions.
#' Default is 0.1 
#' 
#' @param testing An argument for testing purposes which will cause each tumour to be printed as the pliody and mutaitonal GDs are resolved which can aid debugging. 
#' Default is FALSE. 
#' 
#' @param track  An argument which will cause messages to be printed about what the function is doing and a progress bar if a very large amount of data is inputted
#' and the user wishes to track progress. Default is FALSE.
#'
#' @param bootstrap_discovery Whether to assess the bootstrap stability of each 'discover_' threshold call by resampling
#' mutations (with replacement) within each tumour/sample/cluster \code{n_boot} times, and to filter the returned calls on
#' that stability (see @return below for what changes when this is TRUE). Default is FALSE.
#'
#' @param n_boot The number of bootstrap resamples to use when \code{bootstrap_discovery = TRUE}. Default is 500.
#'
#' @param seed The random seed set before bootstrapping, for reproducibility. Default is 1.
#'
#' @param stability_threshold The stability (fraction of bootstrap resamples reproducing a call) a tumour/sample/cluster
#' must clear to be considered stable, used both for the Wilson confidence interval check and as the null value in the
#' binomial significance test. Default is 0.7.
#'
#' @param stability_alpha The significance level used for the Wilson confidence interval on stability, and (when
#' \code{stability_use_binom_test = TRUE}) the q-value cutoff below which a call is considered stable. Default is 0.05.
#'
#' @param stability_use_binom_test Whether stability filtering uses the Benjamini-Hochberg adjusted q-value from a one-sided binomial
#' test (q_value < stability_alpha) rather than the Wilson lower confidence bound (ci_low > stability_threshold).
#' Default is TRUE.
#'
#' @param stability_adjust_method The p-value adjustment method (passed to \code{stats::p.adjust}) used to compute
#' q_value when \code{stability_use_binom_test = TRUE}. Default is "BH".
#'
#'
#' @return A list is returned with three objects which each describe the genome
#' doubling events over all the tumours that were inputted:
#' * GDs_per_tumour: Description of genome doubling events for each tumour (one row per tumour)
#'     * First_GD: Is there a first GD event anywhere in the tumour (ie an event that would modify ploidy from 2 to 4)
#'     * Second_GD: Is there a second GD event anywhere in the tumour (ie an event that would modify ploidy from 3-4 to 6-8)
#'     * num_first_gd: number of first genome doublings across the tumour (clonal or subclonal)
#'     * num_second_gd: number of second genome doublings across the tumour (clonal or subclonal)
#'     * num_clonal_gd: number of clonal genome doublings across the tumour (first or second)
#'     * num_subclonal_gd: number of subclonal genome doublings across the tumour (first or second)
#'     * First_GD_homogen: Whether all regions have any first GD event (from ploidy 2 to 4) 
#'     * Second_GD_homogen: Whether all regions have any second GD event (from ploidy 3-4 to 6-8) 
#'     * GD_status_homogen: Whether all regions have had the same number of GD events (ie will have ~ the same ploidy)
#'     * GD_statuses: A common separated list of all the GD states present in the tumour (from 0, 1 or 2)
#'     * frac_0_gd_regions: The fraction of regions with no GD event
#'     * frac_1_gd_regions:  The fraction of regions with a first GD event (from ploidy 2 to 4)
#'     * frac_2_gd_regions:   The fraction of regions with a second GD event (from ploidy 3-4 to 6-8)
#'     
#' * GDs_per_region: Description of genome doubling events for each region (one row per region)
#'     * gd_clusters: A common separated list of the subclonal clusters in a given region which were identified with 
#'       enough mutations at copy number 2 that some of the mutations probably occurred before a subclonal GD event
#'     * gd_events: A common separated list of the GD event IDs present in a given region
#'     * First_GD: The GD event ID for the first GD (from ploidy 2 to 4) in a given region or 'No GD' if there was no first GD 
#'     * Second_GD: The GD event ID for the second GD (from ploidy 3-4 to 6-8) in a given region or 'No GD' if there was no second GD 

#' * GDs_events: Description of genome doubling events for each tumour (one rwo per event)
#'     * GD_event: The event ID unique within each tumour
#'     * GD_event_id: The event ID concatonated with the tumour id hence is unique accross a cohort
#'     * tumour_id: A id for the tumour
#'     * is_clonal: Whether the event is clonal (present in all regions) or subclonal (pesent in a subset of regions)
#'     * num_regions_gd_event: Number of regions in which the GD is present
#'     * total_regions: Total number of regions in the tumour which the GD event is present
#'     * is_subclonal_mutation_supported: Whether there is a doubled subclonal mutation cluster supporting a subclonal
#'       GD (TRUE) or if this subclonal GD is inferred only from the ploidy (FALSE)
#'     * clusters: If there are supporting subclonal mutation clusters which are they (common seperated list, otherwise NA if none)
#'
#' * mut_counts: The mutation counts per tumour/sample/cluster for the clusters which
#'   support a called subclonal GD event (one row per tumour/sample/cluster which is in the `clusters` field of an event
#'   in GDs_events with `is_subclonal_mutation_supported == TRUE`)
#'
#' * mut_counts_not_called: The same mutation counts as mut_counts, but for every
#'   tumour/sample/cluster combination where the mutation-based test did NOT call a subclonal GD (`is_subcl_gd == FALSE`).
#'   This is intended for diagnosing why a subclonal GD was not called in a given sample/cluster (e.g. too few mutations,
#'   perc_cn2/num_cn2 below threshold etc.) and includes the thresholds used for the run:
#'     * is_subcl_gd: Whether this tumour/sample/cluster was called as showing evidence of a subclonal GD
#'     * is_subcl_gd_any_region: Whether this cluster was called as showing evidence of a subclonal GD in any region/sample
#'       of the tumour (in which case the lower 'check_' thresholds were applied instead of the 'discover_' thresholds)
#'     * discover_mut_cpn_2_threshold, check_mut_cpn_2_threshold, discover_num_muts_threshold, discover_frac_2_cpn_muts_threshold,
#'       discover_num_2_cpn_muts_threshold, check_frac_2_cpn_muts_threshold, check_num_2_cpn_muts_threshold: The threshold values
#'       used for this run (see arguments above)
#'
#' When \code{bootstrap_discovery = TRUE}, calls are filtered through a two-pass procedure before any of the tables
#' above are built, and the function returns extra elements:
#' * GDs_per_tumour, GDs_per_region, GDs_events, mut_counts, mut_counts_not_called: as described above, but built from
#'   a filtered call set:
#'     1. Pass 1 ('anchors'): a tumour/sample/cluster must pass the 'discover_' thresholds AND be stable under
#'        bootstrapping (see \code{stability_use_binom_test}) to be called on its own evidence.
#'     2. Pass 2 ('check', second chance): once a cluster has at least one stable anchor somewhere in the tumour,
#'        every other sample of that same cluster is given a second chance at the lower 'check_' thresholds, without
#'        needing to pass the stability test itself. If nowhere in the tumour is the cluster a significant, stable
#'        subclonal GD, no anchor exists and no other sample of that cluster gets checked. Anchors keep their own
#'        already-validated call (it is not re-evaluated against the check thresholds).
#'   A subclonal GD event can still be called from a single region's stable evidence even when other regions of the
#'   same cluster fail discovery and have no other anchor to borrow from.
#' * unstable_calls_GDs_per_tumour, unstable_calls_GDs_per_region, unstable_calls_GDs_events, unstable_calls_mut_counts,
#'   unstable_calls_mut_counts_not_called: the same tables, but built from whatever the un-filtered discover+check
#'   pipeline would have called, minus the filtered calls above - ie the calls (including any only reached via the
#'   'check' phase with no stable anchor backing the cluster) that are dropped specifically because of stability
#'   filtering - useful to check if you agree with rejected calls.
#' * bootstrap_discovery: One row per tumour/sample/cluster which passed the 'discover_' thresholds on the original
#'   (non-bootstrapped) data, with:
#'     * k, B: number of bootstrap resamples (out of B) in which the call was reproduced
#'     * stability, ci_low, ci_high: k/B and its Wilson confidence interval (level set by \code{stability_alpha})
#'     * keep_call: whether ci_low > stability_threshold
#'     * p_value, q_value: one-sided binomial test (H0: true stability <= stability_threshold) and its
#'       \code{stability_adjust_method}-adjusted q-value (only when \code{stability_use_binom_test = TRUE})
#'     * is_stable: the rule actually used to build the filtered tables above - \code{q_value < stability_alpha}
#'       when \code{stability_use_binom_test = TRUE}, otherwise \code{keep_call}
#' * bootstrap_discovery_settings: the bootstrap/stability parameters used for this run
#'
#' @author
#' 
#' Alexander M Frankell, Francis Crick institute, University College London, \email{alexander.frankell@@crick.ac.uk}
#' Bootstrapping added by Piotr Pawlik, https://github.com/wlippa
#' @examples 
#' # Run on example data (loaded with package)
#' output <- detect_par_gd( example_data )
#' 
#' @export
detect_par_gd <- function( input, discover_mut_cpn_2_threshold = 1.5, check_mut_cpn_2_threshold = NULL,
                           discover_num_muts_threshold = 10,
                           discover_frac_2_cpn_muts_threshold = 0.25, check_frac_2_cpn_muts_threshold = 0.1,
                           discover_num_2_cpn_muts_threshold = 5, check_num_2_cpn_muts_threshold = 3,
                           testing = FALSE, track = FALSE,
                           bootstrap_discovery = FALSE, n_boot = 500, seed = 1,
                           stability_threshold = 0.7, stability_alpha = 0.05,
                           stability_use_binom_test = TRUE, stability_adjust_method = "BH"){

  if (is.null(check_mut_cpn_2_threshold)) {
    check_mut_cpn_2_threshold <- discover_mut_cpn_2_threshold
  }
  
  # get the oriingal class (in case not a data table - revert back at the end) 
  orig_class <- class(input)
  
  # make sure its a data.table for processing
  input <- data.table::as.data.table( input )

  # Remove rows with missing values in mut_cpn, num_gds, or MajCN:
  na_rows <- input[, is.na(mut_cpn) | is.na(num_gds) | is.na(MajCN)]
  if( any(na_rows) ){
    na_samples <- input[ na_rows, unique(sample_id) ]
    warning( sprintf('Removed %d mutations with missing mut_cpn, num_gds or MajCN (samples: %s)',
                     sum(na_rows), paste(na_samples, collapse = ', ')), call. = FALSE )
    input <- input[ !na_rows ]
  }

  input_raw <- data.table::copy(input)
  if( input[, all(num_gds == 0)] ){
    message( 'No GD samples inputted')
    return(NULL)
  }
  
  
  
  if(track) message( 'Detecting evidence of Subclonal GDs from mutations' )
  
  ## for each clone in each region estimate whether at least some of the mutations
  ## were might have been present before a GD event (at mutCPN 2)
  input[, `:=`(num_muts = sum(round(MajCN) == 2^num_gds),
               perc_cn2 = sum(mut_cpn > discover_mut_cpn_2_threshold & round(MajCN) == 2^num_gds) / sum(round(MajCN) == 2^num_gds),
               num_cn2 = sum(mut_cpn > discover_mut_cpn_2_threshold & round(MajCN) == 2^num_gds),
               num_cn2_all = sum(mut_cpn > discover_mut_cpn_2_threshold),
               perc_cn2_all_check = sum(mut_cpn > check_mut_cpn_2_threshold) / .N,
               num_cn2_all_check = sum(mut_cpn > check_mut_cpn_2_threshold),
               num_muts_present = sum(mut_cpn > 0 & round(MajCN) == 2^num_gds),
               num_total_present = sum(mut_cpn > 0) ),
        by = .(tumour_id, cluster_id, sample_id) ]
  input[, is_subcl_gd_discovery := num_muts > discover_num_muts_threshold &
          perc_cn2 > discover_frac_2_cpn_muts_threshold &
          num_cn2 > discover_num_2_cpn_muts_threshold ]
  input[, is_subcl_gd := is_subcl_gd_discovery]

  # set lower threshold ('check_' preflex parameters) if cluster is already doubled in a different region
  # Also remove need for X number of mutations at least at MajCN^num_gds - this aviods calling subclonal GD
  # where in fact different clones have different numbers of mutations at MajCN == numgd^2 (some not enough power for detection)
  input[, is_subcl_gd_any_region := any(is_subcl_gd), by = .(tumour_id, cluster_id)]
  input[ (is_subcl_gd_any_region), is_subcl_gd := perc_cn2_all_check > check_frac_2_cpn_muts_threshold &
                                                  num_cn2_all_check > check_num_2_cpn_muts_threshold]
  
  # thresholds shared by every call to .pgdd_build_outputs below
  threshold_args <- list(
    discover_mut_cpn_2_threshold = discover_mut_cpn_2_threshold,
    check_mut_cpn_2_threshold = check_mut_cpn_2_threshold,
    discover_num_muts_threshold = discover_num_muts_threshold,
    discover_frac_2_cpn_muts_threshold = discover_frac_2_cpn_muts_threshold,
    check_frac_2_cpn_muts_threshold = check_frac_2_cpn_muts_threshold,
    discover_num_2_cpn_muts_threshold = discover_num_2_cpn_muts_threshold,
    check_num_2_cpn_muts_threshold = check_num_2_cpn_muts_threshold,
    testing = testing, track = track
  )

  # is_subcl_gd as computed above already includes the 'check' phase (lower thresholds applied
  # to clusters doubled elsewhere in the tumour). This is the pipeline used as-is when no
  # bootstrapping is requested, and is also the reference ('what would be called with no
  # stability filtering') used below to work out which of those check-phase-rescued calls
  # are 'unstable'.
  full_is_subcl_gd <- input$is_subcl_gd

  if (!isTRUE(bootstrap_discovery)) {

    full_result <- do.call(.pgdd_build_outputs, c(list(input = input, is_subcl_gd_vec = full_is_subcl_gd), threshold_args))
    return( full_result )

  }

  if (!is.numeric(n_boot) || length(n_boot) != 1 || is.na(n_boot) || n_boot < 1) {
    stop("n_boot must be a positive integer.")
  }
  n_boot <- as.integer(n_boot)
  set.seed(seed)

  baseline_discovery <- unique(input[is_subcl_gd_discovery == TRUE,
                                     .(tumour_id, sample_id, cluster_id, num_muts, perc_cn2, num_cn2)])
  if (nrow(baseline_discovery) > 0) {
    baseline_discovery[, `:=`(k = 0L, B = n_boot)]

    for (b in seq_len(n_boot)) {
      boot_input <- input_raw[, .SD[sample.int(.N, .N, replace = TRUE)],
                              by = .(tumour_id, sample_id, cluster_id)]

      boot_disc <- .pgdd_compute_discovery_calls(
        dt = boot_input,
        discover_mut_cpn_2_threshold = discover_mut_cpn_2_threshold,
        discover_num_muts_threshold = discover_num_muts_threshold,
        discover_frac_2_cpn_muts_threshold = discover_frac_2_cpn_muts_threshold,
        discover_num_2_cpn_muts_threshold = discover_num_2_cpn_muts_threshold
      )

      baseline_discovery[boot_disc,
                         on = .(tumour_id, sample_id, cluster_id),
                         k := k + as.integer(i.is_subcl_gd_discovery)]
    }

    ci <- .pgdd_wilson_ci_vec(baseline_discovery$k, baseline_discovery$B, alpha = stability_alpha)
    baseline_discovery[, `:=`(
      stability = k / B,
      ci_low = ci$low,
      ci_high = ci$high
    )]
    baseline_discovery[, keep_call := ci_low > stability_threshold]

    # 'is_stable' is the actual rule used to filter calls below: the q-value from the one-sided
    # binomial test (H0: true stability <= stability_threshold) when stability_use_binom_test is
    # TRUE, or the Wilson lower CI bound ('keep_call') otherwise.
    if (isTRUE(stability_use_binom_test)) {
      baseline_discovery[, p_value := stats::pbinom(k - 1L, size = B, prob = stability_threshold, lower.tail = FALSE)]
      baseline_discovery[, q_value := stats::p.adjust(p_value, method = stability_adjust_method)]
      baseline_discovery[, is_stable := q_value < stability_alpha]
    } else {
      baseline_discovery[, is_stable := keep_call]
    }

    bootstrap_discovery_summary <- baseline_discovery[]
  } else {
    bootstrap_discovery_summary <- data.table::data.table(
      tumour_id = character(), sample_id = character(), cluster_id = character(),
      num_muts = numeric(), perc_cn2 = numeric(), num_cn2 = numeric(),
      k = integer(), B = integer(), stability = numeric(),
      ci_low = numeric(), ci_high = numeric(), keep_call = logical(), is_stable = logical()
    )
  }

  # Pass 1 ('anchors'): a cluster/sample must pass the 'discover_' thresholds AND be stable under
  # bootstrapping to be trusted on its own evidence.
  input[, is_stable := FALSE]
  if (nrow(bootstrap_discovery_summary) > 0) {
    input[bootstrap_discovery_summary, is_stable := i.is_stable,
          on = .(tumour_id, sample_id, cluster_id)]
  }
  input[, anchor_is_subcl_gd := is_subcl_gd_discovery & is_stable]

  # Pass 2 ('check', second chance): once a cluster has at least one stable anchor somewhere in
  # the tumour, every OTHER sample of that same cluster gets a second chance at the lower
  # 'check_' thresholds, without needing to pass the stability test itself - if nowhere in the
  # tumour is the cluster a significant, stable subclonal GD, there is no anchor to justify
  # looking again with lower thresholds. Anchors keep their own already-validated call.
  input[, cluster_has_anchor := any(anchor_is_subcl_gd), by = .(tumour_id, cluster_id)]
  stable_is_subcl_gd <- input[, fifelse(
    anchor_is_subcl_gd, TRUE,
    fifelse(cluster_has_anchor,
            perc_cn2_all_check > check_frac_2_cpn_muts_threshold & num_cn2_all_check > check_num_2_cpn_muts_threshold,
            FALSE) )]
  input[, c('is_stable', 'anchor_is_subcl_gd', 'cluster_has_anchor') := NULL]

  # Unstable calls: whatever the full discover+check pipeline would have called, minus the
  # filtered calls above (anchors plus their check-phase-rescued cluster-mates) - ie the calls
  # (including any only reached via the 'check' phase with no stable anchor backing the cluster)
  # that are lost specifically because of stability filtering.
  unstable_is_subcl_gd <- full_is_subcl_gd & !stable_is_subcl_gd

  stable_result <- do.call(.pgdd_build_outputs, c(list(input = input, is_subcl_gd_vec = stable_is_subcl_gd), threshold_args))
  unstable_result <- do.call(.pgdd_build_outputs, c(list(input = input, is_subcl_gd_vec = unstable_is_subcl_gd), threshold_args))

  # output as list: the default fields carry only stability-filtered ('stable') calls;
  # the calls dropped by stability filtering are kept alongside under unstable_calls_*
  output <- list(GDs_per_tumour = stable_result$GDs_per_tumour,
                 GDs_per_region = stable_result$GDs_per_region,
                 GDs_events = stable_result$GDs_events,
                 mut_counts = stable_result$mut_counts,
                 mut_counts_not_called = stable_result$mut_counts_not_called,
                 unstable_calls_GDs_per_tumour = unstable_result$GDs_per_tumour,
                 unstable_calls_GDs_per_region = unstable_result$GDs_per_region,
                 unstable_calls_GDs_events = unstable_result$GDs_events,
                 unstable_calls_mut_counts = unstable_result$mut_counts,
                 unstable_calls_mut_counts_not_called = unstable_result$mut_counts_not_called,
                 bootstrap_discovery = bootstrap_discovery_summary,
                 bootstrap_discovery_settings = data.table::data.table(
                   n_boot = n_boot,
                   seed = seed,
                   stability_threshold = stability_threshold,
                   stability_alpha = stability_alpha,
                   stability_use_binom_test = stability_use_binom_test,
                   stability_adjust_method = stability_adjust_method,
                   bootstrap_unit = "tumour_id x sample_id x cluster_id",
                   phase = "discovery_only",
                   check_phase_requires_stable_anchor_in_cluster = TRUE
                 ) )

  return( output )

}



seperate_gd_events <- function( tumour_gd_clusters ){
  
  # If no GDs in tumour can return the table back with extra cols
  if(tumour_gd_clusters[, all(num_gds == 0)]){
    tumour_gd_clusters[, gd_events := '' ]
    return( list(tumour_gd_clusters, NULL) )
  }
  
  # get a vector of all doubled clusters in this tumour across all regions
  clusters <- tumour_gd_clusters[, unlist(tstrsplit(gd_clusters, split = ','))]
  clusters <- unique(clusters[ !is.na(clusters) ])
  
  # There may be evidenc of GDs from the pliody but no GD clusters (GD very soon after the MRCA)
  # IN which case call subclonal GDs as normal presuming they are all the same event. If not then do the
  # below
  if( !is.null(clusters) ){
    
    # Do a group clusters based on which are GD'd in the same region (theoretically should represent each event)
    cluster_matrix <- as.data.table( do.call(rbind, lapply(clusters, function(cluster) grepl(cluster, tumour_gd_clusters$gd_clusters)) ))
    cluster_matrix[, num_regions := apply(cluster_matrix, 1, function(x) sum(as.numeric(x))) ]
    colnames( cluster_matrix )[ 1:(ncol(cluster_matrix)-1) ] <-  tumour_gd_clusters$sample_id
    cluster_matrix[, regions_present := apply(cluster_matrix[, 1:(ncol(cluster_matrix)-1)], 1, function(row) paste(row, collapse = ','))]
    cluster_matrix[, cluster := clusters]
    cluster_matrix[, GD_event := .GRP, by = regions_present ]
    
    # Give each of these groupings an 'GD event ID' and sort them by the number of
    # regions which each GD event is in (used for nesting later)
    mut_gd_matrix <- cluster_matrix[, clusters := paste(cluster, collapse = ','), by = GD_event ]
    mut_gd_matrix <- mut_gd_matrix[, `:=`(cluster = NULL, regions_present = NULL) ]
    mut_gd_matrix <- unique( mut_gd_matrix )
    mut_gd_matrix <- mut_gd_matrix[ order(num_regions, decreasing = TRUE) ]
    mut_gd_matrix <- mut_gd_matrix[, GD_event := factor(GD_event, levels = GD_event[order(num_regions,decreasing = TRUE)]) ]
    
    # Need this event table with info per region later
    mut_gd_matrix_save <- mut_gd_matrix
    
    # Get just a matrix of regions and GD events 
    mut_gd_matrix_long <- melt(mut_gd_matrix, id.vars = c('num_regions', 'GD_event', 'clusters'), variable.name = 'region' )
    mut_gd_matrix <- dcast(mut_gd_matrix_long, region ~ GD_event)
    colnames(mut_gd_matrix)[2:ncol(mut_gd_matrix)] <- paste0('GD_event_',colnames(mut_gd_matrix)[2:ncol(mut_gd_matrix)])
    
    # Now compare the number of events to the number of GDs that is indicating by the ploidy
    mut_gd_matrix_org <- copy(mut_gd_matrix)
    mut_gd_matrix[, num_mut_gds := apply(mut_gd_matrix[, 2:ncol(mut_gd_matrix)], 1, function(x) sum(as.numeric(x)))]
    mut_gd_matrix[, num_gds := tumour_gd_clusters[ match(mut_gd_matrix$region, tumour_gd_clusters$sample_id), num_gds] ]
    mut_gd_matrix[, extra_gd := num_gds > num_mut_gds ]
    
    #### ensure never more gds than the ploidy would indicate ####
    mut_gd_matrix[ , extra_mut_gds := num_mut_gds - num_gds ]
    event_cols <-  names(mut_gd_matrix)[ grepl('GD_event', names(mut_gd_matrix)) ]
    
    
    # First deal with regions where more GD events have been identified by the above code that are indicated 
    # by the ploidy. This might because a GD cluster has been missed in one region where it is really present
    # or a large number of amplifications have occured at mutated loci that mimic the affect of GD
    extra_mut_gd_regs <- mut_gd_matrix[ extra_mut_gds > 0, region]
    
    # For these regions remove the GD events so that the number matches that indicated by the pliody  
    if( length(extra_mut_gd_regs) > 0 ){
      
      # Deal with each region wit an extra GD separately
      for(reg in extra_mut_gd_regs){
        
        # how many extra GDs are there for this region
        extra_gds <- mut_gd_matrix[ region == reg, extra_mut_gds]
        
        # What are the mutationally detected GDs in this region?
        events_in_reg <- event_cols[  as.logical( mut_gd_matrix[ region == reg, ..event_cols ]) ]
        num_events <- length(events_in_reg)
        
        # Remove the final events in the table (those in the fewest regions as previously sorted)
        events_to_remove <- events_in_reg[ (num_events - extra_gds + 1):num_events ]
        mut_gd_matrix[ region == reg, (events_to_remove) := FALSE ]
        
      }
      mut_gd_matrix[, extra_mut_gds := NULL ]
      
      # Recalculate the number of GDs 
      mut_gd_matrix[, num_mut_gds := apply(mut_gd_matrix[, 2:(ncol(mut_gd_matrix)-3)], 1, function(x) sum(as.numeric(x)))]
      
      # note whether fewer GDs than expected given the pliod for next section
      mut_gd_matrix[, extra_gd := num_gds > num_mut_gds ]
      
      # remove 'event's that are no longer in any regions
      absent_events <- event_cols[ apply(mut_gd_matrix[ , ..event_cols ],2, function(col) all(!col)) ]
      if( length(absent_events) > 0) set(mut_gd_matrix, j=absent_events, value = NULL )
      event_cols <- event_cols[ !event_cols %in% absent_events ]
    }
    
  } else {
    
    mut_gd_matrix <- data.table( region = tumour_gd_clusters$sample_id,
                                 num_mut_gds = 0,
                                 num_gds = tumour_gd_clusters$num_gds )
    mut_gd_matrix[, extra_gd := num_gds > num_mut_gds ]
    
    
  }
  
  # Now deal with reigons where the pliod indicates more GD events that can be detected by doubled
  # mutant CPN (common - Subclonal GDs that we detect must occur often soon after MRCA). 
  # Add the GDs indicated by ploidy. These must be different events to those indicated by
  # subclonal clusters (as some of the subclonal mutations must have been present before the doubling)
  
  # Use a repeat loop to remove events until these match (usually only one or no repeats needed)
  rep = 0
  while( mut_gd_matrix[, any(extra_gd) ]){
    
    # If no clusters can just presume all regions with the same pliody GD status were part of the same event
    # Otherwise seperate events using the doubled clusters
    if( !is.null(clusters) ){
      
      # First check if any event partially overlaps with another - this nesting structure is impossible for inpdendant events in a phylogeny
      # If so need to remove some events from some regions to make plausible (rare but could occur
      # if many amplifications mimicked a GD event or a mutational GD event was missed)
      event_cols <-  names(mut_gd_matrix)[ grepl('GD_event', names(mut_gd_matrix)) ]
      is_split <- apply( mut_gd_matrix[, ..event_cols ], 2, function(col) any( any(col[ mut_gd_matrix$extra_gd ]) & any(col[ !mut_gd_matrix$extra_gd ]) ) )
      
      # output the additional GD events once nesting and pliody vs mutational GDs are resolved
      if( any(is_split) ){
        split_gd_cols <- event_cols[is_split]
        splits <- as.data.table(apply( mut_gd_matrix[, ..split_gd_cols ], 2, function(col){
          df <- data.table( split_event = col, extra_gds = mut_gd_matrix$extra_gd )
          df[, unique_states := paste(split_event, extra_gds, sep = ',')]
          return( df[ (extra_gds), grp := .GRP, by = unique_states ][, grp ] )
        } ))
        splits[, unique_states := apply(splits, 1, function(row) paste(row, collapse = ',')) ]
        splits[ , grp := .GRP, by = unique_states ]
        splits[ !mut_gd_matrix$extra_gd, grp := NA ]
        splits <- splits[, grp]
      } else {
        
        splits <- rep(1, nrow(mut_gd_matrix) )
        splits[ !mut_gd_matrix$extra_gd ] <- NA
        
      }
      
    } else {
      
      splits <- rep(1, nrow(mut_gd_matrix) )
      splits[ !mut_gd_matrix$extra_gd ] <- NA
      
    } 
    
    # add these additional GD events onto the results table so it is coherent
    for( gd_event in splits[ !is.na(splits) ] ){
      mut_gd_matrix[, (paste0('GD_event_10', gd_event + ((rep)*10))) := splits == gd_event & !is.na(splits) ]
    }
    event_cols <-  names(mut_gd_matrix)[ grepl('GD_event', names(mut_gd_matrix)) ]
    gds <- mut_gd_matrix[, ..event_cols]
    mut_gd_matrix[, num_mut_gds := apply(gds, 1, function(x) sum(as.numeric(x)))]
    mut_gd_matrix[, extra_gd := num_gds > num_mut_gds ]
    
    #repeat if some regions remain unresolved
    rep = rep + 1
  }
  
  # rename all the events in there nesting order (in earliest event = GD1 and later events are GD'>1')
  event_cols <-  names(mut_gd_matrix)[ grepl('GD_event', names(mut_gd_matrix)) ]
  gds <- mut_gd_matrix[, ..event_cols ]
  setcolorder(gds, order(colSums(gds), decreasing = TRUE))
  GD_names <- paste0('GD', 1:ncol(gds))
  orig_names <- colnames(gds)
  names(orig_names) <- GD_names
  setnames(gds, paste0('GD', 1:ncol(gds)))
  
  ## Now create an output where each GD is a row in a table (like mutation table for GDs) with clonality etc indicated ##
  ## First indicate for each region exactly which GD events occured by our estimation
  gds_mat <- do.call(cbind, lapply(1:ncol(gds), function(coli){
    out <- rep('', nrow(gds))
    out[ gds[[coli]] ] <- GD_names[ coli ] 
    return(out)
  }) )
  gd_events <- apply(gds_mat, 1, function(row) paste(row, collapse = ','))
  gd_events <- gsub(',,|,,,|,,,,', ',', gd_events)
  gd_events <- gsub(',,|,,,|,,,,', ',', gd_events)
  gd_events <- gsub('^,|,$', '', gd_events)
  
  # Overlay this onto the region level output
  tumour_gd_clusters[, gd_events := gd_events ]
  
  # now return the table where each GD event is a row
  cols <- c('region', event_cols)
  gd_events <- dcast(melt(mut_gd_matrix[, ..cols ], id.vars = "region", variable.name = 'GD_event'), GD_event ~ region)
  gd_events[, num_regions_gd_event := rowSums(gd_events[, 2:ncol(gd_events)])]
  gd_events[, total_regions := (ncol(gd_events)-2) ]
  gd_events[, samples := sapply(1:nrow(gd_events), function(rowi) paste( names(gd_events)[2:(ncol(gd_events)-2)][as.logical(gd_events[rowi, 2:(ncol(gd_events)-2)])], collapse = ',' ) ) ]
  gd_events[, is_clonal := apply( gd_events[, 2:(ncol(gd_events)-3) ], 1, function(row) all(row) ) ]
  gd_events[, tumour_id := tumour_gd_clusters[, unique(tumour_id) ]]
  
  if( !is.null(clusters) ){
    
    mut_gd_matrix_save[, GD_event := paste0('GD_event_', GD_event)]
    gd_events[, is_subclonal_mutation_supported := GD_event %in% mut_gd_matrix_save$GD_event]
    gd_events[, clusters := mut_gd_matrix_save[ match(gd_events$GD_event, mut_gd_matrix_save$GD_event), clusters ]]
    
  } else {
    
    gd_events[, is_subclonal_mutation_supported := FALSE ]
    gd_events[, clusters := NA ]
    
  }
  gd_events[, GD_event :=  names(orig_names)[ match(GD_event, orig_names)]]
  gd_events[, GD_event_id := paste(gsub('^.{1}_', '', tumour_id), GD_event, sep = '_')]
  
  # return the two different data formats - per region and per GD event
  return( list(tumour_gd_clusters, gd_events[, .(GD_event, GD_event_id, tumour_id, is_clonal, num_regions_gd_event,
                                                 total_regions, is_subclonal_mutation_supported, clusters, samples)] ) )

}

# Takes the per-mutation `input` table (with is_subcl_gd_discovery, is_subcl_gd_any_region and the
# per-cluster diagnostic columns already computed) plus a chosen `is_subcl_gd_vec` (one boolean per
# row of `input`, aligned by row order) and runs the event-resolution / per-tumour and per-region
# summary logic against that specific calling of is_subcl_gd. Used to build the standard output
# tables once for the default (discover+check) pipeline, or twice - once for stability-filtered
# ('stable') calls and once for the calls lost to that filtering ('unstable') - when bootstrapping
# is requested.
.pgdd_build_outputs <- function(input, is_subcl_gd_vec,
                                discover_mut_cpn_2_threshold, check_mut_cpn_2_threshold,
                                discover_num_muts_threshold, discover_frac_2_cpn_muts_threshold,
                                check_frac_2_cpn_muts_threshold,
                                discover_num_2_cpn_muts_threshold, check_num_2_cpn_muts_threshold,
                                testing = FALSE, track = FALSE){

  input <- data.table::copy(input)
  input[, is_subcl_gd := is_subcl_gd_vec]

  # order by most numerous gd clusters (in most samples) - used later to resolve
  input[, num_regions := sum(is_subcl_gd), by = .(tumour_id, cluster_id)]

  # overlay the doubled clusteres for each region
  input[, gd_clusters := paste(unique(cluster_id[(is_subcl_gd & !is_clonal_cluster)]), collapse = ','),
        by = .(sample_id, tumour_id)]

  #### Now need to work out what is the simplest explanation of events to lead to these clusters
  #### being genome doubled ####
  # Reduce table to per region and cluster GDs / pliody GDs
  input_subcl_clusters <- unique( input[, .(gd_clusters, num_gds), by = .(sample_id, tumour_id)] )

  if(track) message( 'Resolving with pliody and nesting structure for each tumour' )

  tumours <- input_subcl_clusters[, unique(tumour_id)]
  if(track) pb <- utils::txtProgressBar( min = 0, max = length(tumours), style = 3, width =  30 )
  mut_gds_both <- lapply( tumours, function(tumour){

    if(testing) print(tumour)
    if(track) utils::setTxtProgressBar( pb, which( tumours == tumour ) )
    seperate_gd_events( tumour_gd_clusters = input_subcl_clusters[ tumour_id == tumour ] )

  }  )
  mut_gds_seperated <-  rbindlist( lapply(mut_gds_both, function(x) x[[1]]) )
  mut_gds_events <- rbindlist( lapply(mut_gds_both, function(x) x[[2]]) )
  mut_gds_events <- mut_gds_events[ order(GD_event_id) ]

  # Clean up the NAs
  mut_gds_events[ is.na(clusters), clusters := NA ]

  # add mutation count:
  # include all mutation count inputs used for threshold-based subclonal GD calling.
  relevant_columns = c('tumour_id', 'sample_id', 'cluster_id', 'num_muts', 'num_cn2',
                       'num_cn2_all', 'num_muts_present', 'perc_cn2',
                       'num_cn2_all_check', 'perc_cn2_all_check',
                       'num_total_present', 'is_subcl_gd_discovery', 'is_subcl_gd', 'is_subcl_gd_any_region')
  mut_counts_all =  unique(input[, ..relevant_columns])
  mut_counts_all[, cluster_id := as.character(cluster_id) ]

  # Add explicit threshold diagnostics per row so users can see:
  # value tested, threshold used, and whether each threshold was cleared.
  mut_counts_all[, `:=`(
    discover_mut_cpn_2_threshold = discover_mut_cpn_2_threshold,
    check_mut_cpn_2_threshold = check_mut_cpn_2_threshold,

    discover_num_muts_value = num_muts,
    discover_num_muts_threshold = discover_num_muts_threshold,
    discover_num_muts_pass = num_muts > discover_num_muts_threshold,

    discover_frac_2_cpn_muts_value = perc_cn2,
    discover_frac_2_cpn_muts_threshold = discover_frac_2_cpn_muts_threshold,
    discover_frac_2_cpn_muts_pass = perc_cn2 > discover_frac_2_cpn_muts_threshold,

    discover_num_2_cpn_muts_value = num_cn2,
    discover_num_2_cpn_muts_threshold = discover_num_2_cpn_muts_threshold,
    discover_num_2_cpn_muts_pass = num_cn2 > discover_num_2_cpn_muts_threshold,

    check_frac_2_cpn_muts_value = perc_cn2_all_check,
    check_frac_2_cpn_muts_threshold = check_frac_2_cpn_muts_threshold,
    check_frac_2_cpn_muts_pass = perc_cn2_all_check > check_frac_2_cpn_muts_threshold,

    check_num_2_cpn_muts_value = num_cn2_all_check,
    check_num_2_cpn_muts_threshold = check_num_2_cpn_muts_threshold,
    check_num_2_cpn_muts_pass = num_cn2_all_check > check_num_2_cpn_muts_threshold
  )]

  mut_counts_all[, `:=`(
    discover_pass_all = discover_num_muts_pass &
      discover_frac_2_cpn_muts_pass &
      discover_num_2_cpn_muts_pass,
    check_pass_all = check_frac_2_cpn_muts_pass &
      check_num_2_cpn_muts_pass
  )]
  if (all(is.na(mut_gds_events$clusters))) {
    clusters_supporting_wgd = c()
    } else {
  clusters_supporting_wgd = mut_gds_events[ is_subclonal_mutation_supported == TRUE, unlist(strsplit(clusters, split = ',')) ]
    }
  mut_counts = mut_counts_all[ cluster_id %in% clusters_supporting_wgd ]

  # Diagnostic table: every tumour/sample/cluster combination for which the mutation-based
  # test did NOT call a subclonal GD (is_subcl_gd == FALSE), together with the underlying
  # counts and the thresholds used in this run. This is intended to help work out why a
  # subclonal GD was not called for a given sample/cluster (e.g. too few mutations,
  # perc_cn2/num_cn2 below threshold etc.)
  mut_counts_not_called <- mut_counts_all[ is_subcl_gd == FALSE ]

  # Summarise per tumour
  mut_gds_seperated[, First_GD := tstrsplit(gd_events, split = ',')[[1]]]
  if( mut_gds_seperated[, any( grepl(',', gd_events) )]){
    mut_gds_seperated[, Second_GD := tstrsplit(gd_events, split = ',')[[2]]]
  } else {
    mut_gds_seperated[, Second_GD := as.character(NA) ]
  }

  mut_gds_seperated[ is.na(First_GD), First_GD := 'No GD' ]
  mut_gds_seperated[ is.na(Second_GD), Second_GD := 'No GD' ]

  mut_gds_tumour <- mut_gds_seperated[, .(First_GD = ifelse( any(!First_GD == 'No GD'), ifelse(length(unique(First_GD)) > 1, 'Subclonal', 'Clonal'), 'No GD'),
                                          Second_GD = ifelse( any(!Second_GD == 'No GD'), ifelse(length(unique(Second_GD)) > 1, 'Subclonal', 'Clonal'), 'No GD'),
                                          num_first_gd = length(unique(First_GD[ !First_GD == 'No GD' ])),
                                          num_second_gd = length(unique(Second_GD[ !Second_GD == 'No GD' ])),
                                          First_GD_homogen = all(First_GD == unique(First_GD)[1]),
                                          Second_GD_homogen = all(Second_GD == unique(Second_GD)[1]),
                                          GD_status_homogen = all(num_gds == unique(num_gds)[1]),
                                          GD_statuses = paste(unique(num_gds)[ order(unique(num_gds)) ], collapse = ','),
                                          frac_0_gd_regions = sum(num_gds == 0)/.N,
                                          frac_1_gd_regions = sum(num_gds == 1)/.N,
                                          frac_2_gd_regions = sum(num_gds == 2)/.N),
                                      by = tumour_id ]

  # over the number of clonal and subclonal gds calculated from the per event table
  mut_gds_events_tumour <- mut_gds_events[, .( num_clonal_gds = sum(is_clonal == TRUE),
                                              num_subclonal_gds = sum(is_clonal == FALSE)),
                                          by = tumour_id ]
  mut_gds_tumour <- merge(mut_gds_tumour, mut_gds_events_tumour, by = 'tumour_id', all.x = TRUE)
  mut_gds_tumour[ is.na(num_clonal_gds), num_clonal_gds := 0 ]
  mut_gds_tumour[ is.na(num_subclonal_gds), num_subclonal_gds := 0 ]

  list(GDs_per_tumour = mut_gds_tumour,
      GDs_per_region = mut_gds_seperated,
      GDs_events = mut_gds_events,
      mut_counts = mut_counts,
      mut_counts_not_called = mut_counts_not_called)
}

.pgdd_compute_discovery_calls <- function(dt,
                                          discover_mut_cpn_2_threshold,
                                          discover_num_muts_threshold,
                                          discover_frac_2_cpn_muts_threshold,
                                          discover_num_2_cpn_muts_threshold) {
  dt <- data.table::as.data.table(dt)

  out <- dt[, .(
    num_muts = sum(round(MajCN) == 2^num_gds),
    perc_cn2 = {
      denom <- sum(round(MajCN) == 2^num_gds)
      if (denom == 0) 0 else sum(mut_cpn > discover_mut_cpn_2_threshold & round(MajCN) == 2^num_gds) / denom
    },
    num_cn2 = sum(mut_cpn > discover_mut_cpn_2_threshold & round(MajCN) == 2^num_gds)
  ), by = .(tumour_id, cluster_id, sample_id)]

  out[, is_subcl_gd_discovery := num_muts > discover_num_muts_threshold &
        perc_cn2 > discover_frac_2_cpn_muts_threshold &
        num_cn2 > discover_num_2_cpn_muts_threshold]
  out
}

.pgdd_wilson_ci_vec <- function(k, n, alpha = 0.05) {
  z <- stats::qnorm(1 - alpha / 2)
  phat <- k / n
  denom <- 1 + (z^2 / n)
  centre <- (phat + (z^2 / (2 * n))) / denom
  half <- (z / denom) * sqrt((phat * (1 - phat) / n) + (z^2 / (4 * n^2)))
  list(low = pmax(0, centre - half), high = pmin(1, centre + half))
}

#############
###  END  ###
#############


discover_mut_cpn_2_threshold = 1.5; check_mut_cpn_2_threshold = 1.25; discover_num_muts_threshold = 10;
discover_frac_2_cpn_muts_threshold = 0.25; check_frac_2_cpn_muts_threshold = 0.1;
discover_num_2_cpn_muts_threshold = 5; check_num_2_cpn_muts_threshold = 3;
testing = FALSE; track = FALSE
