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
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character"),
    make_option(c("-r", "--rep"), type = "integer", default = 1, 
                help = "Replicate number.  Special case 0 means use base (for real data).", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile
rep <- opt$rep

# go where we want data outputs to be
setwd( '../data' )
setwd( dir_out )
if ( rep > 0 ) {
    dir_rep <- paste0( 'rep-', rep )
    setwd( dir_rep )
}

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
    softlink( 'std_mor', name, 'grm.id' ) # these files exist but are bad, overwrite!
    softlink( 'std_mor', name, 'grm.N.bin' )
}

# There are all of these kinship matrices to consider
# (true and *_lim available on simulation only)
# - `kinship_true`: true kinship matrix of simulation
# - `kinship_popkin_rom`: unbiased (but noisy) estimate from genotypes
# - `kinship_popkin_mor`: MOR version of popkin, possibly biased but performs better on GWAS
# - `kinship_std_mor`: biased (and noisy) MOR version of "standard" estimate from genotypes (same as GEMMAs)
# - `kinship_std_rom`: biased (and noisy) ROM version of "standard" estimate from genotypes
# - `kinship_std_rom_lim`: limit of biased ROM version of "standard" estimator
# - `kinship_wg_mor`: MOR version of WG (based on popkin MOR)
# - `kinship_wg_rom`: biased (and noisy) WG estimate from genotypes
# - `kinship_wg_rom_lim`: limit of biased WG estimator

# behavior depends on the presence of a true kinship matrix, which tells us if this is a simulation or a real dataset.
# present => sim; absent => real.
is_sim <- file.exists( 'kinship/true.grm.bin' )

# in a real dataset, this directory hasn't been created yet, so do it now!
if ( !is_sim )
    dir.create( 'kinship' )

# standard estimate, calculated with GCTA
# (only case not calculated in R, need gcta64 binary)
# do first as other estimators link aux files to these
gcta_grm( name, name_out = 'kinship/std_mor' )
delete_files_log( 'kinship/std_mor' ) # unneeded log file

# rest of work occurs in this subdirectory
setwd( 'kinship' )

# biased "standard" kinship estimate, ROM only
write_grm_softlinks( 'std_rom', kinship_std( X ) )
#write_grm_softlinks( 'std_mor', kinship_std( X, mean_of_ratios = TRUE ) ) # same as GCTA's

# popkin estimates, without labels
kinship_popkin_rom <- popkin( X )
write_grm_softlinks( 'popkin_rom', kinship_popkin_rom )
kinship_popkin_mor <- popkin( X, mean_of_ratios = TRUE )
write_grm_softlinks( 'popkin_mor', kinship_popkin_mor )

# WG estimate
write_grm_softlinks( 'wg_rom', kinship_wg_limit( kinship_popkin_rom ) )
# get this new MOR version for free from popkin!
write_grm_softlinks( 'wg_mor', kinship_wg_limit( kinship_popkin_mor ) )

# limits are avaiable for simulations only (true kinship must be known)
if ( is_sim ) {
    # true kinship (as given by simulation; already written but need to load to calculate limits of biased estimates)
    kinship_true <- read_grm( 'true' )$kinship / 2
    # (re)link these optional files to STD's so they're identical (and save space)
    softlink( 'std_mor', 'true', 'grm.id' )
    softlink( 'std_mor', 'true', 'grm.N.bin' )

    # limit of biased "standard" estimator
    write_grm_softlinks( 'std_rom_lim', kinship_std_limit( kinship_true ) )

    # WG limit
    write_grm_softlinks( 'wg_rom_lim', kinship_wg_limit( kinship_true ) )
}

