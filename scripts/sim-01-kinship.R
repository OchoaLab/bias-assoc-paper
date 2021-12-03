# calculate various kinship matrices (estimated and limits)
library(optparse)    # for terminal options
library(BEDMatrix)
library(genio)       # to write files for external software
library(popkin)      # to estimate kinship without bias
library(popkinsuppl) # for PCA's kinship estimator
library(genbin)      # gcta and plink binary wrappers

# a name for temporary BED/etc data, under project dir
name <- 'data'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile

# go where we want data outputs to be
setwd( '../data' )
setwd( dir_out )

# load genotypes for kinship estimates
X <- BEDMatrix( name, simple_names = TRUE )

###############
### KINSHIP ###
###############

# helper functions

# this version overwrites new file if it existed
softlink <- function( name1, name2, ext ) {
    file1 <- paste0( name1, '.', ext )
    file2 <- paste0( name2, '.', ext )
    system3( 'ln', c( '-f', '-s', file1, file2 ) )
}

# automated version for most cases
write_grm_softlinks <- function(name, kinship) {
    # write main GRM files
    write_grm( name, 2 * kinship )
    # create a softlink from true "M" file, so those are the same for all files (limit variation due to that)
    # same with ID table (less big but whatever)
    # softlink saves space!
    softlink( 'gcta', name, 'grm.id' ) # these files exist but are bad, overwrite!
    softlink( 'gcta', name, 'grm.N.bin' )
}

# There are all of these kinship matrices to consider
# - `kinship_true`: true kinship matrix of simulation
# - `kinship_std_rom`: biased (and noisy) ROM version of "standard" estimate from genotypes
# - `kinship_std_rom_lim`: limit of biased ROM version of "standard" estimator
# - `kinship_std_mor`: biased (and noisy) MOR version of "standard" estimate from genotypes (same as GEMMAs)
# - `kinship_popkin`: unbaised (but noisy) estimate from genotypes
# - `kinship_wg`: biased (and noisy) WG estimate from genotypes
# - `kinship_wg_lim`: limit of biased WG estimator
# - `kinship_gcta`: biased (and noisy) GCTA estimate from genotypes
# - `kinship_gcta_lim`: limit of biased GCTA estimator

# behavior depends on the presence of a true kinship matrix, which tells us if this is a simulation or a real dataset.
# present => sim; absent => real.
is_sim <- file.exists( 'kinship/true.grm.bin' )

# in a real dataset, this directory hasn't been created yet, so do it now!
if ( !is_sim )
    dir.create( 'kinship' )

# GCTA estimate
# (only case not calculated in R, need gcta64 binary)
# do first as other estimators link aux files to these
gcta_grm( name, name_out = 'kinship/gcta' )
delete_files_log( 'kinship/gcta' ) # unneeded log file

# rest of work occurs in this subdirectory
setwd( 'kinship' )

# biased "standard" kinship estimate
write_grm_softlinks( 'std_rom', kinship_std( X ) )
write_grm_softlinks( 'std_mor', kinship_std( X, mean_of_ratios = TRUE ) )

# popkin estimate, without labels
kinship_popkin <- popkin( X )
write_grm_softlinks( 'popkin', kinship_popkin )

# WG estimate
write_grm_softlinks( 'wg', kinship_wg_limit( kinship_popkin ) )

# limits are avaiable for simulations only (true kinship must be known)
if ( is_sim ) {
    # true kinship (as given by simulation; already written but need to load to calculate limits of biased estimates)
    kinship_true <- read_grm( 'true' )$kinship / 2
    # (re)link these optional files to GCTA's so they're identical (and save space)
    softlink( 'gcta', 'true', 'grm.id' )
    softlink( 'gcta', 'true', 'grm.N.bin' )

    # limit of biased "standard" estimator
    write_grm_softlinks( 'std_rom_lim', kinship_std_limit( kinship_true ) )

    # WG limit
    write_grm_softlinks( 'wg_lim', kinship_wg_limit( kinship_true ) )

    # GCTA limit
    write_grm_softlinks( 'gcta_lim', kinship_gcta_limit( kinship_true ) )
}

