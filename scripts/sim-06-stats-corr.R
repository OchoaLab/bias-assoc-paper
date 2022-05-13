library(optparse) # for terminal options
library(readr)    # to read tables
library(popkin)   # to plot
library(ochoalabtools) # for nice PDF
library(tibble)

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character"),
    make_option(c("-r", "--rep"), type = "integer", default = 1, 
                help = "Replicate number", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile
rep <- opt$rep

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cc' )

# load pre-existing data
setwd( '../data/' )
setwd( dir_out )

# if this directory exists here, this is real data
is_sim <- !file.exists( 'kinship' )

dir_rep <- paste0( 'rep-', rep )
setwd( dir_rep )

# in real data the oracle methods (truth and biased limits) are missing, exclude from `kinship_methods`
if ( !is_sim )
    kinship_methods <- kinship_methods[ kinship_methods$type != 'ROM lim.', ]

# get correct count of methods now
n_kinship <- nrow( kinship_methods )

# let's plot data in the order of the `kinship_methods` table
# also, plot PCA first, LMM second
method_codes <- c(
    paste0( 'pca_', kinship_methods$code ),
    paste0( 'lmm_', kinship_methods$code )
)
# same but human-readable
labs_short <- rep.int( kinship_methods$short, 2 )
labs_type <- rep.int( kinship_methods$type, 2 )
labs_model <- c(
    rep.int( 'PCA', n_kinship ),
    rep.int( 'LMM', n_kinship )
)

# load tibbles
pvals <- read_tsv( 'pvals.txt.gz', col_types = cols( ) )
betas <- read_tsv( 'betas.txt.gz', col_types = cols( ) )

plot_cor <- function( data, name ) {
    # reorder columns of data to be a manually-selected order
    data <- data[ method_codes ]
    # replace original codes with nice names
    colnames( data ) <- labs_short
    # compute correlation matrix
    cor_data <- cor( data, use = 'pairwise.complete.obs' )
    
    # get nice max width for a journal
    dim <- fig_width()
    fig_start(
        name,
        width = dim,
        height = dim * 0.84
    )
    plot_popkin(
        cor_data,
        labs = cbind( labs_model, labs_type ),
        labs_cex = c(1, 0.7),
        labs_line = c(6, 4),
        labs_even = TRUE,
        names = TRUE,
        names_cex = 0.7,
        mar = 7,
        ylab = 'Association Model, Kinship Estimate',
        ylab_adj = 0.75,
        leg_title = expression(bold(paste("Pearson Correlation ", (rho)))),
        leg_width = 0.15
    )
    fig_end()

    # return in same order, etc
    return( cor_data )
}

# save figure in lower level
setwd( '..' )

cor_pvals <- plot_cor( pvals, 'pvals_cor' )
# plot of betas fail in TGP (min correlations are much too high), just skip
if ( is_sim )
    cor_betas <- plot_cor( betas, 'betas_cor' )

# we'll see which of these subsets are of interest in real data (were manually picked for sim data only)
if ( is_sim ) {
    # pick out some values of particular importance
    # p-values only
    range_subset <- function( indexes ) {
        # take subset
        cor_pvals <- cor_pvals[ indexes, indexes ]
        # return row, to grow a tibble with data
        tibble(
            subset = toString( method_codes[ indexes ] ),
            min = min( cor_pvals ),
            max = max( cor_pvals )
        )
    }

    data <- NULL
    # these should be all PCA methods except for True and Popkin
    data <- rbind( data, range_subset( c( 2:4, 6:9 ) ) )
    # include Popkin too
    data <- rbind( data, range_subset( 2:9 ) )
    # include True too
    data <- rbind( data, range_subset( 1:9 ) )
    # now LMM subsets
    # limits only except GCTA
    data <- rbind( data, range_subset( 10:12 ) )
    # all limits only
    data <- rbind( data, range_subset( 10:13 ) )
    # estimates ROM: popkin-wg-std
    data <- rbind( data, range_subset( 14:16 ) )
    # estimates MOR: std-gcta
    data <- rbind( data, range_subset( 17:18 ) )
    # all estimators only
    data <- rbind( data, range_subset( 14:18 ) )
    # all LMM
    data <- rbind( data, range_subset( 10:18 ) )

    # save report
    write_tsv( data, 'pvals_cor_report.txt'  )
}
