# just calculate V matrices, save them!
library(optparse)    # for terminal options
library(readr)
library(genio) # for read_grm
library(simtrait) # for cov_trait

# wrapper around this setup
read_kin <- function( name_method )
    read_grm( paste0( 'kinship/', name_method ) )$kinship / 2

# these use global `dir_herit`
# this one has no scaling and other such bells and whistles
write_V <- function( V, name_method )
    write_grm( paste0( dir_herit, 'V/', name_method ), V )

calc_V_single <- function( data, rep, kinships, name_method ) {
    # determine if output already exists, do not recalculate!
    file_out <- paste0( dir_herit, 'V/', name_method, '.grm.bin' )
    if ( file.exists( file_out ) ) {
        message( 'Output exists, skipping: ', file_out )
        return()
    }
    
    # extract sigmas for this method from the data
    # first reverse-engineer the type and short name code from the full name
    if ( name_method == 'true' ) {
        type <- 'rom_lim'
        short <- 'tr'
    } else {
        type <- sub( '^(popkin|std|wg)_', '', name_method )
        short <- sub( paste0( '_', type, '$'), '', name_method )
        if ( short == 'popkin' ) short <- 'tr'
    }
    # this combination yields a single row!
    data <- data[ data$rep == rep & data$type == type, ]
    # now get the columns we want
    sigma_sq_g <- data[[ paste0( 'g_', short ) ]]
    sigma_sq_e <- data[[ paste0( 'e_', short ) ]]

    # load matrix from rep, unless provided (for real data, as they are shared across reps)
    if ( is.null( kinships ) )
        kinship <- read_kin( name_method )
    else
        kinship <- kinships[[ name_method ]]
    
    # calculate V matrix
    sigma_sq_p <- sigma_sq_g + sigma_sq_e
    herit <- sigma_sq_g / sigma_sq_p
    # most accurate fit of this matrix, with scale included
    V <- cov_trait( kinship, herit, sigma_sq_p )
    
    # write this matrix to an output file
    write_V( V, name_method )
}

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character"),
    make_option("--n_rep", type = "integer", default = NA, 
                help = "Total number of replicates", metavar = "int"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile
n_rep <- opt$n_rep
herit <- opt$herit

if ( is.na( n_rep ) )
    stop( 'Option `--n_rep` is required!' )

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cccc' )

# go where we want data outputs to be
setwd( '../data' )
setwd( dir_out )

# if this directory exists here, this is real data
is_sim <- !file.exists( 'kinship' )

# to get back down easily
dir_base <- getwd()

# include additional level if heritability is non-default
dir_herit <- '' # so stuff works for default case too
if ( herit != 0.8 ) {
    dir_herit <- paste0( 'h-', herit, '/' )
    setwd( dir_herit )
}

# load sigma!
data <- read_tsv( 'reml.txt.gz', col_types = 'icdddddddd' )

# for non-default herit, have to get back to base!
# (this works for default herit too)
setwd( dir_base )

# NULL means load kinship matrices from each rep as we go
kinships <- NULL
kinship_names <- kinship_methods$code
# since real only has one set of each kinship estimate, load them upfront
if ( !is_sim ) {
    # no limits in this case
    kinship_names <- kinship_methods$code[ kinship_methods$type != 'ROM lim.' ]
    # load all upfront
    kinships <- lapply( kinship_names, read_kin )
    # these names don't get added automatically, do it here!
    names( kinships ) <- kinship_names
}

# get back down if needed
setwd( dir_base )

# navigate all replicates
for ( rep in 1 : n_rep ) {
    message( 'rep: ', rep )
    dir_rep <- paste0( 'rep-', rep )
    setwd( dir_rep )
    # make output dir for this rep
    dir_out_rep <- paste0( dir_herit, 'V' )
    if ( !dir.exists( dir_out_rep ) )
        dir.create( dir_out_rep )
    
    # process each type, calculating and saving its V
    lapply( kinship_names, function( name ) calc_V_single( data, rep, kinships, name ) )
    
    setwd( '..' )
}

