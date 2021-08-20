
#default variables
covariates = c("ID","entity","gender")
write_intermediate=T
#test modes: 1. fixed genomic bin, 2. not implemented:genes
test_mode = 1	
minentityposcasespct = 0.05
maxentityposcasespct = 1.01
thresholdhypmutation=100000  #~10mut/Mb of 3,234.83 Mb
min_testing_bin_vars = 6
framesize_pruned = 30
frame_pruned_min_nvar = 6
corrVariantFactor = 1.5