# ParallelGDDetect
A package to detect distinct genome doubling events which have occurred in different samples from the same tumour.

A function is provided which takes as input the number of genome doublings (GDs) for each region as estimated from the genome
wide copy number and mutation copy numbers in each region and devolves which GDs across samples are part of
the same event and which present distinct events as indicated by the doubling of subclonal mutations which
will occur when the subclonal mutation arises before a given subclonal GD event. We recommend plotting the mutation
copy numbers alongside the allele specific copy number across the genome to verify these calls as well as carefully
checking the ploidy solutions determined for each sample. Default thresholds have been set for best performance
on TRACER NSCLC exome sequencing data. Particularly for whole genome sequencing data, thresholds
numbers may need modification.

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
