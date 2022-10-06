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
                help = "Total number of replicates", metavar = "int"),
    make_option("--noWG", action = "store_true", default = FALSE, 
                help = "Create plot version that excludes WG, for some presentations")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile
n_rep <- opt$n_rep
noWG <- opt$noWG
if ( is.na( n_rep ) )
    stop( 'Option `--n_rep` is required!' )

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cccc' )

# go where the data is
setwd( '../data/' )
setwd( dir_out )

# if this directory exists here, this is real data
is_sim <- !file.exists( 'kinship' )

# load precalculated data
data <- read_tsv( 'eval.txt.gz', col_types = 'cid' )

# in real data the oracle methods (truth and biased limits) are missing, exclude from `kinship_methods`
if ( !is_sim )
    kinship_methods <- kinship_methods[ kinship_methods$type != 'ROM lim.', ]

# exclude WG here!
if ( noWG )
    kinship_methods <- kinship_methods[ grep( '^wg_', kinship_methods$code, invert = TRUE ), ]

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

    # for noWG make output path different
    name_out <- name
    dims <- fig_scale( 1.5 ) # w/h
    if ( noWG ) {
        name_out <- paste0( name_out, '-noWG' )
        # don't want these to be full width, so approach this way instead
        width <- fig_width() / 2
        height <- width
        dims <- c( width, height )
    }
    
    # now make plot of data
    fig_start(
        name_out,
        width = dims[1],
        height = dims[2],
        mar_b = 8
    )
    boxplot(
        data_list,
        names = NA, # individual labels will be plotted with rest
        xaxt = 'n',
        ylab = if ( name == 'rmsd' ) lab_rmsd else lab_auc
    )
    mtext( 'Association Model, Kinship Estimate', side = 1, line = 7 )
    # mark zero and band area in this case
    if ( name == 'rmsd' ) {
        abline( h = 0, lty = 2, col = 'gray' )
        abline( h = 0.01, lty = 2, col = 'gray' )
        abline( h = -0.01, lty = 2, col = 'gray' )
    }
    
    # reuse popkin-style labeling!  Great for hierarchical/factorial setup
    # though popkin has defaults for these, the raw function doesn't have defaults!
    popkin:::print_labels_multi(
                 labs = cbind( labs_short, labs_type, labs_model ),
                 labs_cex = c(0.7, 0.7, 1),
                 labs_las = c(2, 0, 0),
                 labs_line = c(0.5, 4, 6),
                 labs_lwd = 1, # default
                 labs_sep = c(FALSE, TRUE, TRUE),
                 labs_even = c(FALSE, TRUE, TRUE),
                 labs_ticks = FALSE, # default
                 labs_text = TRUE, # default
                 labs_col = 'black', # default
                 # these align with barplot-specific changes
                 xc_ind = 1 : length( method_codes ),
                 doMat = FALSE
             )
    
    fig_end()
}

plot_measure( 'auc' )
plot_measure( 'rmsd' )

