# this script creates a small figure with fewer data, specifically for grant proposals
# intended for TGP only, but won't hardcode for now

library(optparse) # for terminal options
library(readr)    # to read tables
library(ochoalabtools) # for nice PDF
library(popkin)

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
# directly filter data we want to focus on here
codes_small <- c('popkin_rom', 'popkin_mor', 'std_rom', 'std_mor')
kinship_methods <- kinship_methods[ kinship_methods$code %in% codes_small, ]
# edit names some more (space is precious)
kinship_methods$nice <- sub( ' est.', '', kinship_methods$nice )

# go where the data is
setwd( '../data/' )
setwd( dir_out )

# load precalculated data
data <- read_tsv( 'eval.txt.gz', col_types = 'cid' )

# get correct count of methods now
n_kinship <- nrow( kinship_methods )

# let's plot data in the order of the `kinship_methods` table
# LMM only for this case!
method_codes <- paste0( 'lmm_', kinship_methods$code )
# same but human-readable
labs_short <- kinship_methods$short
labs_type <- kinship_methods$type

# plotting labels
lab_auc <- expression( bold( AUC[PR] ) )

# reorganize data for boxplots
# will appear in desired order (from kinship_methods)
data_list <- lapply( method_codes, function( x ) {
    # subset tibble to data from this kinship method only
    auc_x <- data[ data$kinship == x, ]
    # validate number of replicates
    stopifnot( nrow( auc_x ) == n_rep )
    stopifnot( all( auc_x$rep %in% 1 : n_rep ) )
    # just keep column of interest, only values of interest
    return( auc_x$auc )
})
# here we'll go with normal names under boxplots
names( data_list ) <- kinship_methods$nice

# now make plot of data
fig_start(
    'auc-small',
    width = 2.5, # way smaller than full fig, for grant
    height = 3.5,
    mar_b = 5.5
)
# to control labels underneath plot
par_orig <- par( las = 3, cex.axis = 0.7 )
boxplot(
    data_list,
    ylab = lab_auc
)
# restore `las` and anything else that might have gotten messed up
par( par_orig )
mtext( 'Kinship estimate', side = 1, line = 4.5 )

fig_end()
