# ParallelGDDetect
A package to detect distinct genome doubling events which have occurred in different samples from the same tumour.

Function is provided which takes as input the number of genome doublings (GDs) for each region as estimated from the genome
wide copy number and mutation copy numbers in each region and devolves which GDs across samples are part of
the same event and which present distinct events as indicated by the doubling of subclonal mutations which
will occur when the subclonal mutation arises before a given subclonal GD event. We recommend plotting the mutation
copy numbers alongside the allele specific copy number across the genome to verify these calls as well as carefully
checking the ploidy solutions determined for each sample. Default thresholds have been set for best performance
on TRACER NSCLC exome sequencing data. Particularly for whole genome sequencing data, thresholds
numbers may need modification.
