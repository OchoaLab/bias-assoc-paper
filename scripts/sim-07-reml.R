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

data_reml <- NULL

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

    # return at least these
    data <- c( sigma_sq_g, sigma_sq_e )

    if ( factors ) {
        # load kinship
        kinship <- read_grm( name_grm )$kinship / 2
    
        # also calculate expected factors relating the last two to the first
        # the factor is "c" from the paper, and multiplies the first to get the others
        c_st <- 1 - mean( kinship )
        c_wg <- 1 - mean( kinship[ lower.tri( kinship ) ] )
        # append to return vector
        data <- c( data, c_st, c_wg )
    }

    return( data )
}

reml_tab <- function( data_tr, data_st, data_wg, rep, type ) {
    # organize values into tibble for easy output
    data_reml <<- bind_rows(
        data_reml,
        tibble(
            rep = rep,
            type = type,
            g_tr = data_tr[1],
            e_tr = data_tr[2],
            g_std = data_st[1],
            e_std = data_st[2],
            g_wg = data_wg[1],
            e_wg = data_wg[2],
            c_std = data_tr[3],
            c_wg = data_tr[4]
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
        reml_tab( data_tr, data_st, data_wg, rep, 'rom_lim' )
    }

    data_tr <- gcta_reml_sigmas( 'popkin_rom', factors = TRUE )
    data_st <- gcta_reml_sigmas( 'std_rom' )
    data_wg <- gcta_reml_sigmas( 'wg_rom' )
    reml_tab( data_tr, data_st, data_wg, rep, 'rom' )

    data_tr <- gcta_reml_sigmas( 'popkin_mor', factors = TRUE )
    data_st <- gcta_reml_sigmas( 'std_mor' )
    data_wg <- gcta_reml_sigmas( 'wg_mor' )
    reml_tab( data_tr, data_st, data_wg, rep, 'mor' )

    # go back down for simulations (only case where it is needed)
    if ( is_sim ) 
        setwd( '..' )
}

# save data frame!
write_tsv( data_reml, 'reml.txt.gz' )
