library(optparse) # for terminal options
library(readr)    # to read tables
library(ochoalabtools) # for nice PDF
library(dplyr)    # for bind_rows

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
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cc' )

# go where the data is
setwd( '../data/' )
setwd( dir_out )

# if this directory exists here, this is real data
is_sim <- !file.exists( 'kinship' )

# tibble to grow
data <- NULL

# load pre-calculated AUCs and SRMSDs
for ( rep in 1 : n_rep ) {
    file_rep <- paste0( 'rep-', rep, '/eval.txt.gz' )
    data_rep <- read_tsv( file_rep, col_types = 'cid' )
    data <- bind_rows( data, data_rep )
}

if ( !is_sim ) {
    # in real data the oracle methods (truth and biased limits) are missing
    # harmonize `kinship_methods` to match
    # "true" is only one without "_lim" at the end, the rest can be picked up with that regex
    codes_rm <- c('true', grep('_lim$', kinship_methods$code, value = TRUE) )
    # now prune table
    kinship_methods <- kinship_methods[ !(kinship_methods$code %in% codes_rm), ]
    # since all that is left are estimators, could remove the "est." part of the "nice" names
    kinship_methods$nice <- sub( ' est.', '', kinship_methods$nice )
}
# get correct count of methods now
n_kinship <- nrow( kinship_methods )


# let's plot data in the order of the `kinship_methods` table
# also, plot PCA first, LMM second
method_codes <- c(
    paste0( 'pca_', kinship_methods$code ),
    paste0( 'lmm_', kinship_methods$code )
)
# compute AUC for each case
n_methods <- length( method_codes ) # = 2 * n_kinship

# plotting labels
lab_rmsd <- expression( bold( SRMSD[p] ) )
lab_auc <- expression( bold( AUC[PR] ) )

# shared plotting for AUC and RMSD
# uses lots of globals
plot_measure <- function( name = 'auc' ) {
    # reorganize data for boxplots
    # will appear in desired order (from kinship_methods)
    data_list <- lapply( method_codes, function( x ) {
        # subset tibble to data from this kinship method only
        auc_x <- data[ data$kinship == x, ]
        # validate number of replicates
        stopifnot( nrow( auc_x ) == n_rep )
        stopifnot( all( auc_x$rep %in% 1 : n_rep ) )
        # just keep column of interest, only values of interest
        return( auc_x[[ name ]] )
    })

    # replace codes with nice names we want on plot
    # NOTE: no PCA/LMM marks (will be added to plot separately)
    names( data_list ) <- rep.int( kinship_methods$nice, 2 )

    # now make plot of data
    dims <- fig_scale( if ( is_sim ) 1.5 else 2 ) # w/h
    fig_start(
        name,
        width = dims[1],
        height = dims[2],
        mar_l = if ( is_sim ) 11 else 10
    )
    boxplot(
        rev( data_list ),
        horizontal = TRUE,
        las = 1,
        xlab = if ( name == 'rmsd' ) lab_rmsd else lab_auc
    )
    # mark zero and band area in this case
    if ( name == 'rmsd' ) {
        abline( v = 0, lty = 2, col = 'gray' )
        abline( v = 0.01, lty = 2, col = 'gray' )
        abline( v = -0.01, lty = 2, col = 'gray' )
    }
        
    # add separating lines
    x_line0 <- if ( is_sim ) {
                   if ( name == 'rmsd' ) -0.083 else 0.07
               } else {
                   if ( name == 'rmsd' ) -0.025 else -0.005
               }
    x_line_txt <- x_line0 * if ( name == 'rmsd' ) 1.04 else if ( is_sim ) 0.85 else 2
    x_line <- c(x_line0, x_line0)
    y_line1 <- c(1, n_kinship)
    y_line2 <- c(1, n_kinship) + n_kinship
    lines( x_line, y_line1, xpd = NA )
    lines( x_line, y_line2, xpd = NA )
    # and labels
    text( x_line_txt, mean( y_line2 ), 'PCA', xpd = NA, srt = 90 )
    text( x_line_txt, mean( y_line1 ), 'LMM', xpd = NA, srt = 90 )
    fig_end()
}

plot_measure( 'auc' )
plot_measure( 'rmsd' )

