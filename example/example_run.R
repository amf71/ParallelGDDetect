###### script to test that the function is still working properly after any changes ######

#should change to function that automatically gets to git parent dir
setwd("/Volumes/proj-tracerx-lung/tctProjects/frankella/repos/my_R_packages/ParrallelGDDetect/")

# load required libraries
library(data.table)
library(dplyr)

# source the latest code #
source('R/parrallel_gd_caller.R')

# load the example data #

example_files <- system("ls data", intern = TRUE)

for( file in example_files ) load( paste0("data/", file) )

# example_data
# tumour_id chromosome  position ref alt cluster_id is_clonal_cluster sample_id mut_cpn num_gds
# 1:  CRUK0196      chr19  39964012   C   A          1              TRUE  SU_T1-R1    9.98       2
# 2:  CRUK0196      chr19  39964012   C   A          1              TRUE  SU_T1-R3   15.73       2
# 3:  CRUK0196      chr19  39964012   C   A          1              TRUE  SU_T1-R4   13.01       2
# 4:  CRUK0196       chr3 160083937   G   T          1              TRUE  SU_T1-R1    9.09       2
# 5:  CRUK0196       chr3 160083937   G   T          1              TRUE  SU_T1-R3    6.00       2
# ---
#   48609:  CRUK0495       chr6 155571770   C   T         27             FALSE  SU_T1-R4    0.00       1
# 48610:  CRUK0495       chr6 155571770   C   T         27             FALSE  SU_T1-R6    0.00       1
# 48611:  CRUK0495       chr6 155571770   C   T         27             FALSE  SU_T2-R1    0.00       1
# 48612:  CRUK0495       chr6 155571770   C   T         27             FALSE  SU_T2-R2    0.24       1
# 48613:  CRUK0495       chr6 155571770   C   T         27             FALSE  SU_T3-R2    0.00       1

## run the examples and output to test folder ##

output <- detect_par_gd( input = example_data )

# lapply(output, function(x) head(x,8))
# $GDs_per_tumour
# tumour_id  First_GD Second_GD num_first_gd num_second_gd First_GD_homogen Second_GD_homogen GD_status_homogen GD_statuses
# 1:  CRUK0196    Clonal Subclonal            1             1             TRUE             FALSE             FALSE         1,2
# 2:  CRUK0380    Clonal     No GD            1             0             TRUE              TRUE              TRUE           1
# 3:  CRUK0495 Subclonal Subclonal            3             1            FALSE             FALSE             FALSE         1,2
# 4:  CRUK0516    Clonal     No GD            1             0             TRUE              TRUE              TRUE           1
# 5:  CRUK0565 Subclonal     No GD            1             0            FALSE              TRUE             FALSE         0,1
# 6:  CRUK0681 Subclonal     No GD            2             0            FALSE              TRUE             FALSE         0,1
# 7:  CRUK0714    Clonal     No GD            1             0             TRUE              TRUE              TRUE           1
# frac_0_gd_regions frac_1_gd_regions frac_2_gd_regions num_clonal_gds num_subclonal_gds
# 1:         0.0000000         0.3333333         0.6666667              1                 1
# 2:         0.0000000         1.0000000         0.0000000              1                 0
# 3:         0.0000000         0.7500000         0.2500000              0                 4
# 4:         0.0000000         1.0000000         0.0000000              1                 0
# 5:         0.2500000         0.7500000         0.0000000              0                 1
# 6:         0.7142857         0.2857143         0.0000000              0                 2
# 7:         0.0000000         1.0000000         0.0000000              1                 0
# 
# $GDs_per_region
# sample_id tumour_id gd_clusters num_gds gd_events First_GD Second_GD
# 1:  SU_T1-R1  CRUK0495                   1       GD1      GD1     No GD
# 2:  SU_T1-R2  CRUK0495                   1       GD1      GD1     No GD
# 3:  SU_T1-R3  CRUK0495                   1       GD1      GD1     No GD
# 4:  SU_T1-R4  CRUK0495                   1       GD1      GD1     No GD
# 5:  SU_T1-R6  CRUK0495           2       1       GD4      GD4     No GD
# 6:  SU_T2-R1  CRUK0495         5,8       2   GD2,GD3      GD2       GD3
# 7:  SU_T2-R2  CRUK0495         5,8       2   GD2,GD3      GD2       GD3
# 8:  SU_T3-R2  CRUK0495           5       1       GD2      GD2     No GD
# 
# $GDs_events
# GD_event  GD_event_id tumour_id is_clonal num_regions_gd_event total_regions is_subclonal_mutation_supported clusters
# 1:      GD1 CRUK0196_GD1  CRUK0196      TRUE                    3             3                           FALSE     <NA>
#   2:      GD2 CRUK0196_GD2  CRUK0196     FALSE                    2             3                           FALSE     <NA>
#   3:      GD1 CRUK0380_GD1  CRUK0380      TRUE                    8             8                           FALSE     <NA>
#   4:      GD1 CRUK0495_GD1  CRUK0495     FALSE                    4             8                           FALSE     <NA>
#   5:      GD2 CRUK0495_GD2  CRUK0495     FALSE                    3             8                            TRUE        5
# 6:      GD3 CRUK0495_GD3  CRUK0495     FALSE                    2             8                            TRUE        8
# 7:      GD4 CRUK0495_GD4  CRUK0495     FALSE                    1             8                            TRUE        2
# 8:      GD1 CRUK0516_GD1  CRUK0516      TRUE                    8             8                           FALSE     <NA>
#   samples
# 1:                                              SU_T1-R1,SU_T1-R3,SU_T1-R4
# 2:                                                       SU_T1-R1,SU_T1-R3
# 3: SU_T1-R1,SU_T1-R2,SU_T1-R3,SU_T1-R4,SU_T1-R5,SU_T1-R6,SU_T1-R7,SU_T1-R8
# 4:                                     SU_T1-R1,SU_T1-R2,SU_T1-R3,SU_T1-R4
# 5:                                              SU_T2-R1,SU_T2-R2,SU_T3-R2
# 6:                                                       SU_T2-R1,SU_T2-R2
# 7:                                                                SU_T1-R6
# 8: SU_T1-R1,SU_T1-R2,SU_T1-R3,SU_T1-R4,SU_T1-R5,SU_T1-R6,SU_T1-R7,SU_T1-R8


# save
fwrite(output[[1]], 'example/example_output/example_GDs_per_tumour.tsv', sep = '\t')
fwrite(output[[2]], 'example/example_output/example_GDs_per_region.tsv', sep = '\t')
fwrite(output[[3]], 'example/example_output/example_GD_events.tsv', sep = '\t')

#######
# END #
#######























# first simple example #

png( "data-raw/testing/test_outputs/example_1.png" )

cloneMap( tree.mat = tree_example_1, CCF.data = CCFs_example_1, border.thickness = 3 )

invisible( dev.off() )


# second more complex example #

png( "data-raw/testing/test_outputs/example_2.png")

cloneMap( tree_example_2, CCFs_example_2, border.thickness = 3)

invisible( dev.off() )


# second example again using clone_map object #

clone_map_eg <- cloneMap( tree_example_2, CCFs_example_2, output.Clone.map.obj = TRUE, plot.data = FALSE )

png( "data-raw/testing/test_outputs/example_3.png" )

cloneMap( clone_map = clone_map_eg, border.thickness = 3 )

invisible( dev.off() )


# plot both egs as if they were from the same tumours to show how to specify same colours for clones #

layout <- matrix( c( 1, 2,
                     3, 4 ), ncol = 2, byrow = TRUE )

png( "data-raw/testing/test_outputs/example_4.png", width = 800 )

layout( layout,
        heights = c(1, 5),
        widths = c(2, 2))

par( mar = c(1, 1, 1, 1), xpd = NA)

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0.5, 0.5, labels = "Tumour 1: Sample 1", cex = 3, font = 2)

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0.5, 0.5, labels = "Tumour 1: Sample 2", cex = 3, font = 2)

cloneMap( tree_example_1, CCFs_example_1, clone.cols = clone_colours_example, border.thickness = 3 )

cloneMap( tree_example_2, CCFs_example_2, clone.cols = clone_colours_example, border.thickness = 3 )

invisible( dev.off() )



# plot maps of polyclonal data similar to that found in normal tissues

png( "data-raw/testing/test_outputs/example_polyclonal.png", width = 800 )

cloneMap::cloneMap( tree.mat = tree_example_poly, 
                    CCF.data = CCF_example_poly )

invisible( dev.off() )

# plot maps of polyclonal data similar to that found in normal tissues with border

png( "data-raw/testing/test_outputs/example_polyclonal_border.png", width = 800 )

cloneMap::cloneMap( tree.mat = tree_example_poly, 
                    CCF.data = CCF_example_poly,
                    tissue_border = TRUE)

invisible( dev.off() )

# plot maps of polyclonal data similar to that found in normal tissues with border

png( "data-raw/testing/test_outputs/example_polyclonal_spaced.png", width = 800 )

cloneMap::cloneMap( tree.mat = tree_example_poly, 
                    CCF.data = CCF_example_poly,
                    tissue_border = TRUE,
                    space_fraction = 0.7 )

invisible( dev.off() )


### END ###



# for testing min function line by line

# tree.mat = NA; CCF.data = NA; clone_map = NA; output.Clone.map.obj = FALSE;
# plot.data = TRUE; high_qualty_mode = FALSE; track = NA; brewer.palette = "Paired";
# clone.cols = NA; border.colour = "grey20";  border.thickness = 1.5;
# resolution.index = 100;  smoothing.par = 15; repeat.limit = 4; space_fraction = NA;
# tissue_border = FALSE
