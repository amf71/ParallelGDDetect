#====================================================#
#     Function to produce test for Parrallel GDs     #
#    across different samples from the same tumour   #
#====================================================#



# input = example_data
# input = eg_input_all
# 
# mut_cpn_2_threshold = 1.5; discover_num_muts_threshold = 20;
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
#' @param mut_cpn_2_threshold The mutation copy number (mut_cpn) theshold to consider a mutation most likely at 2 copies. Deafult is 1.5. 
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
#' @author 
#' 
#' Alexander M Frankell, Francis Crick institute, University College London, \email{alexander.frankell@@crick.ac.uk}
#' 
#' @examples 
#' # Run on example data (loaded with package)
#' output <- detect_par_gd( example_data )
#' 
#' @export
detect_par_gd <- function( input, mut_cpn_2_threshold = 1.5, discover_num_muts_threshold = 10,
                           discover_frac_2_cpn_muts_threshold = 0.25, check_frac_2_cpn_muts_threshold = 0.1,
                           discover_num_2_cpn_muts_threshold = 5, check_num_2_cpn_muts_threshold = 3, 
                           testing = FALSE, track = FALSE){
  
  # get the oriingal class (in case not a data table - revert back at the end) 
  orig_class <- class(input)
  
  # make sure its a data.table for processing
  input <- data.table::as.data.table( input )
  
  if( input[, all(num_gds == 0)] ){
    message( 'No GD samples inputted')
    return(NULL)
  }
  
  
  
  if(track) message( 'Detecting evidence of Subclonal GDs from mutations' )
  
  ## for each clone in each region estimate whether at least some of the mutations
  ## were might have been present before a GD event (at mutCPN 2)
  input[, `:=`(num_muts = sum(round(MajCN) == 2^num_gds),
               perc_cn2 = sum(mut_cpn > mut_cpn_2_threshold & round(MajCN) == 2^num_gds) / sum(round(MajCN) == 2^num_gds),
               num_cn2 = sum(mut_cpn > mut_cpn_2_threshold & round(MajCN) == 2^num_gds),
               perc_cn2_all = sum(mut_cpn > mut_cpn_2_threshold) / .N,
               num_cn2_all = sum(mut_cpn > mut_cpn_2_threshold) ),
        by = .(tumour_id, cluster_id, sample_id) ]
  input[, is_subcl_gd := num_muts > discover_num_muts_threshold & 
          perc_cn2 > discover_frac_2_cpn_muts_threshold &
          num_cn2 > discover_num_2_cpn_muts_threshold ]
  
  # set lower threshold ('check_' preflex parameters) if cluster is already doubled in a different region
  # Also remove need for X number of mutations at least at MajCN^num_gds - this aviods calling subclonal GD
  # where in fact different clones have different numbers of mutations at MajCN == numgd^2 (some not enough power for detection)
  input[, is_subcl_gd_any_region := any(is_subcl_gd), by = .(tumour_id, cluster_id)]
  input[ (is_subcl_gd_any_region), is_subcl_gd := perc_cn2_all > check_frac_2_cpn_muts_threshold &
                                                  num_cn2_all > check_num_2_cpn_muts_threshold]
  
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
  
  # output as list
  output <- list(GDs_per_tumour = mut_gds_tumour, 
                 GDs_per_region = mut_gds_seperated, 
                 GDs_events = mut_gds_events)
  
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


#############
###  END  ###
#############


mut_cpn_2_threshold = 1.5; discover_num_muts_threshold = 10;
discover_frac_2_cpn_muts_threshold = 0.25; check_frac_2_cpn_muts_threshold = 0.1;
discover_num_2_cpn_muts_threshold = 5; check_num_2_cpn_muts_threshold = 3; 
testing = FALSE; track = FALSE












