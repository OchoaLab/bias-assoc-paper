# this script redoes REML fitting, tests our predictions about sigmas of bias and unbiased kinship matrices
library(optparse)    # for terminal options
library(readr)       # to write output tables
library(genbin)      # gcta and plink binary wrappers
library(genio)
library(simtrait)
library(tibble)
library(dplyr)

# a name for temporary BED/etc data, under project dir
name <- 'data'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character"),
    make_option("--n_rep", type = "integer", default = NA, 
                help = "Total number of replicates", metavar = "int"),
    make_option(c("-t", "--threads"), type = "integer", default = 0, 
                help = "number of threads (default use all cores)", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile
n_rep <- opt$n_rep
threads <- opt$threads
if ( is.na( n_rep ) )
    stop( 'Option `--n_rep` is required!' )

# go where we want data outputs to be
setwd( '../data' )
setwd( dir_out )

# if this directory exists here, this is real data
is_sim <- !file.exists( 'kinship' )

data_err <- NULL

gcta_reml_sigmas <- function( name_method, factors = FALSE ) {
    # setup paths and names for output tables
    name_grm <- paste0( 'kinship/', name_method )

    message( dir_rep, ' ', name_method )
    # REML
    data <- gcta_reml( name, name_grm = name_grm, name_phen = name_phen, threads = threads )
    # cleanup
    delete_files_gcta_hsq( name )
    delete_files_log( name )

    # extract the two sigmas of interest from the existing data
    data <- data$data # want this table only
    sigma_sq_g <- data$Variance[ data$Source == 'V(G)' ]
    sigma_sq_e <- data$Variance[ data$Source == 'V(e)' ]
    # minimal validation
    stopifnot( sigma_sq_g >= 0 )
    stopifnot( sigma_sq_e >= 0 )

    # load kinship for posdef test at least
    kinship <- read_grm( name_grm )$kinship / 2
    
    # also want to answer the question of whether the total covariance for WG ends up being positive definite!
    sigma_sq_p <- sigma_sq_g + sigma_sq_e
    herit <- sigma_sq_g / sigma_sq_p
    # most accurate fit of this matrix, with scale included
    V <- cov_trait( kinship, herit, sigma_sq_p )
    # use min eigenvalue as evidence of non-posdef, just return as is
    emin_k <- min( eigen( kinship, symmetric = TRUE, only.values = TRUE )$values )
    emin_v <- min( eigen( V, symmetric = TRUE, only.values = TRUE )$values )
    
    # return at least these
    data <- c( sigma_sq_g, sigma_sq_e, emin_k, emin_v )

    if ( factors ) {
        # also calculate expected factors relating the last two to the first
        # the factor is "c" from the paper, and multiplies the first to get the others
        c_st <- 1 - mean( kinship )
        c_wg <- 1 - mean( kinship[ lower.tri( kinship ) ] )
        # append to return vector
        data <- c( data, c_st, c_wg )
    }

    return( data )
}

bias_pred_errors <- function( data_tr, data_st, data_wg, rep, type ) {
    # condense the various numbers further, into errors
    
    # for environmental, prediction is that they are all the same
    # for genetic variance, use factor predicted from theory
    # also store all minimum eigenvalues, for later analysis
    # put them directly in a table, smaller than the input!
    data_err <<- bind_rows(
        data_err,
        tibble(
            rep = rep,
            type = type,
            err_e_st = data_st[2] - data_tr[2],
            err_e_wg = data_wg[2] - data_tr[2],
            err_g_st = data_st[1] - data_tr[1] * data_tr[5],
            err_g_wg = data_wg[1] - data_tr[1] * data_tr[6],
            emin_k_tr = data_tr[3],
            emin_v_tr = data_tr[4],
            emin_k_st = data_st[3],
            emin_v_st = data_st[4],
            emin_k_wg = data_wg[3],
            emin_v_wg = data_wg[4]
        )
    )
}

for ( rep in 1 : n_rep ) {
    dir_rep <- paste0( 'rep-', rep )
    # stay in lower level for real data
    if ( is_sim ) 
        setwd( dir_rep )

    # include rep dir for phenotypes if we have real data
    name_phen <- if ( is_sim ) name else paste0( dir_rep, '/', name )

    ############
    ### REML ###
    ############

    # limits are available for simulations only (true kinship must be known)
    if ( is_sim ) {
        data_tr <- gcta_reml_sigmas( 'true', factors = TRUE )
        data_st <- gcta_reml_sigmas( 'std_rom_lim' )
        data_wg <- gcta_reml_sigmas( 'wg_rom_lim' )
        bias_pred_errors( data_tr, data_st, data_wg, rep, 'rom_lim' )
    }

    data_tr <- gcta_reml_sigmas( 'popkin_rom', factors = TRUE )
    data_st <- gcta_reml_sigmas( 'std_rom' )
    data_wg <- gcta_reml_sigmas( 'wg_rom' )
    bias_pred_errors( data_tr, data_st, data_wg, rep, 'rom' )

    data_tr <- gcta_reml_sigmas( 'popkin_mor', factors = TRUE )
    data_st <- gcta_reml_sigmas( 'std_mor' )
    data_wg <- gcta_reml_sigmas( 'wg_mor' )
    bias_pred_errors( data_tr, data_st, data_wg, rep, 'mor' )

    # go back down for simulations (only case where it is needed)
    if ( is_sim ) 
        setwd( '..' )
}

# save data frame!
write_tsv( data_err, 'preds-reml.txt.gz' )
