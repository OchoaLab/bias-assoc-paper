library(optparse) # for terminal options
library(readr)    # to read tables
library(popkin)   # to plot
library(ochoalabtools) # for nice PDF

# constants
# used for legend labels
tolp <- 1e-2
tolb <- 1e-3

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character"),
    make_option("--noWG", action = "store_true", default = FALSE, 
                help = "Create plot version that excludes WG, for some presentations")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile
noWG <- opt$noWG

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cccc' )

# load pre-existing data
setwd( '../data/' )
setwd( dir_out )

# if this directory exists here, this is real data
is_sim <- !file.exists( 'kinship' )

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

# will also save in addition to plot
plot_cor <- function( name, eq = FALSE, tol ) {
    # load precalculated values
    cor_data <- as.matrix( read_tsv( paste0( name, '.txt' ), show_col_types = FALSE ) )

    # for noWG make output path different
    name_out <- name
    # get nice max width for a journal
    width <- fig_width()
    height <- width * 0.84
    leg_width <- 0.15
    if ( noWG ) {
        name_out <- paste0( name_out, '-noWG' )
        width <- width * 2 / 3 # shrink a bit?
        height <- width * 0.75
        leg_width <- 0.25
        # also subset data!
        indexes <- colnames( cor_data ) %in% method_codes
        cor_data <- cor_data[ indexes, indexes ]
    }
    # either way, transfer short names to matrix, they will be used as names
    colnames( cor_data ) <- labs_short
    rownames( cor_data ) <- labs_short
    
    fig_start(
        name_out,
        width = width,
        height = height
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
        leg_title = if (eq) paste0( 'Proportion |diff| < ', tol ) else 'Pearson Correlation',
        leg_width = leg_width
    )
    fig_end()
}

plot_cor( 'pvals_cor' )
plot_cor( 'pvals_eq', eq = TRUE, tol = tolp )
plot_cor( 'betas_eq', eq = TRUE, tol = tolb )
# plot of betas fail in TGP (min correlations are much too high)
if ( is_sim )
    plot_cor( 'betas_cor' )
