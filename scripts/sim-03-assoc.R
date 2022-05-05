# runs association tests (PCA and LMM) iterating over kinship estimates
library(optparse)    # for terminal options
library(readr)       # to write output tables
library(genio)       # to write BED files for external software
library(genbin)      # gcta and plink binary wrappers
library(tibble)

# a name for temporary BED/etc data, under project dir
name <- 'data'
# same name for all cases! (gets overwritten serially, removed in the end)
file_covar <- paste0( name, '.eigenvec' )

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character"),
    make_option(c("-r", "--n_pcs"), type = "integer", default = 2, 
                help = "Number of PCs to use", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile
n_pcs <- opt$n_pcs

# indexes for PCs
indexes <- 1 : n_pcs

# go where we want data outputs to be
setwd( '../data' )
setwd( dir_out )

# behavior depends on the presence of a true kinship matrix, which tells us if this is a simulation or a real dataset.
# present => sim; absent => real.
is_sim <- file.exists( 'kinship/true.grm.bin' )

#############
### ASSOC ###
#############

# get number of loci from plink data, needed to initialize tibbles
m_loci <- count_lines( name, ext = 'bim' )
# empty tibbles to grow
pvals <- tibble( .rows = m_loci )
betas <- tibble( .rows = m_loci )

# automated version for all other cases (only GCTA with default kinship matrix makes more sense differently)
assoc_all <- function( name_method ) {
    # setup paths and names for output tables
    name_grm <- paste0( 'kinship/', name_method )
    name_lmm_method <- paste0( 'lmm_', name_method ) # for LMM outputs
    name_pca_method <- paste0( 'pca_', name_method ) # for PCA outputs

    # LMM
    message( name_lmm_method )
    # association
    data <- gcta_mlma( name, name_grm = name_grm )
    # cleanup
    delete_files_gcta_mlma( name )
    delete_files_log( name )
    # add to data frames
    pvals[[ name_lmm_method ]] <- data$p
    betas[[ name_lmm_method ]] <- data$beta
    
    # PCA
    message( name_pca_method )
    # have to load kinship matrix first
    data <- read_grm( name_grm )
    kinship <- data$kinship / 2
    fam <- data$fam
    eigenvectors <- eigen( kinship )$vectors
    eigenvectors <- eigenvectors[, indexes ] # subset
    # write eigenvectors to file
    write_eigenvec( name, eigenvectors, fam, plink2 = TRUE )
    # set maximum VIF for two special cases where near-collinearity is otherwise a problem for plink with defaults
    # Actual error message: "  Error: Cannot proceed with --glm regression on phenotype 'PHENO1', since variance inflation factor for covariate 'PC1' is too high (VIF_TOO_HIGH). You may want to remove redundant covariates and try again."
    vif <- if ( name_method == 'true' || grepl( 'popkin', name_method ) ) 100 else NA
    # association
    data <- plink_glm( name, file_covar = file_covar, vif = vif )
    # cleanup
    unlink( file_covar )
    delete_files_plink_glm( name )
    delete_files_log( name )
    # add to data frames
    pvals[[ name_pca_method ]] <- data$p
    betas[[ name_pca_method ]] <- data$beta
    
    # make sure globals are edited
    pvals <<- pvals
    betas <<- betas
}

# limits are avaiable for simulations only (true kinship must be known)
if ( is_sim ) {
    assoc_all( 'true' )
    assoc_all( 'std_rom_lim' )
    assoc_all( 'gcta_rom_lim' )
    assoc_all( 'wg_rom_lim' )
}

assoc_all( 'popkin_rom' )
assoc_all( 'popkin_mor' )
assoc_all( 'std_rom' )
assoc_all( 'std_mor' )
assoc_all( 'gcta_mor' )
assoc_all( 'wg_rom' )
assoc_all( 'wg_mor' )

# save data frames!
write_tsv( pvals, 'pvals.txt.gz' )
write_tsv( betas, 'betas.txt.gz' )
