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
