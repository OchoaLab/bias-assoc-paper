# just calculate V matrices, save them!
library(optparse)    # for terminal options
library(readr)
library(genio) # for read_grm
library(simtrait) # for cov_trait
library(tibble)
library(dplyr) # for bind_rows

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character"),
    make_option("--n_rep", type = "integer", default = NA, 
                help = "Total number of replicates", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile
n_rep <- opt$n_rep
if ( is.na( n_rep ) )
    stop( 'Option `--n_rep` is required!' )

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cccc' )

# go where we want data outputs to be
setwd( '../data' )
setwd( dir_out )

# wrapper around this setup
read_kin_process <- function( name_method, rep, is_V = FALSE ) {
    # load matrix
    type <- if ( is_V ) 'V' else 'kinship'
    kinship <- read_grm( paste0( type, '/', name_method ) )$kinship
    # kinship matrices are halved, V's don't need that
    if ( !is_V ) kinship <- kinship / 2

    # eigendecompose! (slowest part)
    evs <- eigen( kinship, symmetric = TRUE, only.values = TRUE )$values
    # get the two stats we want
    # first is minimum eigenvalue
    emin <- min( evs )
    # second is condition number
    evs <- abs( evs )
    kappa <- max( evs ) / min( evs )
    # gather in a tibble row, return
    return( tibble( kinship = name_method, rep = rep, type = type, emin = emin, kappa = kappa ) )
}

# tibble to grow
data <- NULL

# if this directory exists here, this is real data
is_sim <- !file.exists( 'kinship' )

# NULL means load kinship matrices from each rep as we go
kinship_names <- kinship_methods$code
# since real only has one set of each kinship estimate, process them once
if ( !is_sim ) {
    # no limits in this case
    kinship_names <- kinship_methods$code[ kinship_methods$type != 'ROM lim.' ]
    # process all
    for ( kinship_name in kinship_names ) {
        data_rep <- read_kin_process( kinship_name, rep = 0L )
        data <- bind_rows( data, data_rep )
    }
}

# navigate all replicates
for ( rep in 1L : n_rep ) {
    message( 'rep: ', rep )
    dir_rep <- paste0( 'rep-', rep )
    setwd( dir_rep )

    # process all
    for ( kinship_name in kinship_names ) {
        if ( is_sim ) {
            # process kinship matrices in reps for sims only
            data_rep <- read_kin_process( kinship_name, rep = rep )
            data <- bind_rows( data, data_rep )
        }
        # always process V matrices in reps
        data_rep <- read_kin_process( kinship_name, rep = rep, is_V = TRUE )
        data <- bind_rows( data, data_rep )
    }
    
    setwd( '..' )
}

# write table to output
write_tsv( data, 'eigen.txt.gz' )
