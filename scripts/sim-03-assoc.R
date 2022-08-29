# runs association tests (PCA and LMM) iterating over kinship estimates
library(optparse)    # for terminal options
library(readr)       # to write output tables
library(genio)       # to write BED files for external software
library(genbin)      # gcta and plink binary wrappers
library(tibble)

# increase default VIF threshold to this much, for true K and popkin only, admix-fam-sim only
vif_true_popkin_sim <- 20000
vif_true_popkin_real <- NA # 150
vif_default <- NA
max_corr_true_popkin_sim <- 1
max_corr_default <- NA
# max_corr_true_popkin_real is NA (set through if/else logic below)
# --max-corr default 0.999, actual 0.9990483 (rep-24, the bad one), -0.9856413 (rep-1), 0.9282595 for same sim w/ G=1

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
                help = "Replicate number", metavar = "int"),
    make_option(c("-p", "--n_pcs"), type = "integer", default = 2, 
                help = "Number of PCs to use", metavar = "int"),
    make_option(c("-t", "--threads"), type = "integer", default = 0, 
                help = "number of threads (default use all cores)", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile
rep <- opt$rep
n_pcs <- opt$n_pcs
threads <- opt$threads

# indexes for PCs
indexes <- 1 : n_pcs

# go where we want data outputs to be
setwd( '../data' )
setwd( dir_out )

# if this directory exists here, this is real data
is_sim <- !file.exists( 'kinship' )

dir_rep <- paste0( 'rep-', rep )
# stay in lower level for real data
if ( is_sim ) 
    setwd( dir_rep )

# include rep dir for phenotypes if we have real data
name_phen <- if ( is_sim ) name else paste0( dir_rep, '/', name )

# to run in parallel on cluster, make sure outputs (even temporary ones) from different reps do not overlap!
name_out <- paste0( name, '-rep-', rep )
# same name for all cases! (gets overwritten serially, removed in the end)
file_covar <- paste0( name_out, '.eigenvec' )

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
    data <- gcta_mlma( name, name_grm = name_grm, name_phen = name_phen, name_out = name_out, threads = threads )
    # cleanup
    delete_files_gcta_mlma( name_out )
    delete_files_log( name_out )
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
    write_eigenvec( name_out, eigenvectors, fam, plink2 = TRUE )
    # set maximum VIF for two special cases where near-collinearity is otherwise a problem for plink with defaults
    # Actual error message: "  Error: Cannot proceed with --glm regression on phenotype 'PHENO1', since variance inflation factor for covariate 'PC1' is too high (VIF_TOO_HIGH). You may want to remove redundant covariates and try again."
    vif <- if ( name_method == 'true' || grepl( 'popkin', name_method ) ) { if (is_sim) vif_true_popkin_sim else vif_true_popkin_real } else vif_default
    # ditto maximum correlation, observed for true/popkin in admix-fam-sim only
    max_corr <- if ( is_sim && ( name_method == 'true' || grepl( 'popkin', name_method ) ) ) max_corr_true_popkin_sim else max_corr_default
    # association
    data <- plink_glm( name, name_phen = name_phen, file_covar = file_covar, vif = vif, max_corr = max_corr, name_out = name_out, threads = threads )
    # cleanup
    unlink( file_covar )
    delete_files_plink_glm( name_out )
    delete_files_log( name_out )
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
    assoc_all( 'wg_rom_lim' )
}

assoc_all( 'popkin_rom' )
assoc_all( 'popkin_mor' )
assoc_all( 'std_rom' )
assoc_all( 'std_mor' )
assoc_all( 'wg_rom' )
assoc_all( 'wg_mor' )

# save data frames!
# move into rep now for real data (we're already there for sims)
if ( !is_sim )
    setwd( dir_rep )
write_tsv( pvals, 'pvals.txt.gz' )
write_tsv( betas, 'betas.txt.gz' )
