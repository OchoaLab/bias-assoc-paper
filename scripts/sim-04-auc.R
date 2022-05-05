library(optparse) # for terminal options
library(readr)    # to read tables
library(ochoalabtools) # for nice PDF
library(simtrait) # pval_aucpr

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cc' )

# load pre-existing data
setwd( '../data/' )
setwd( dir_out )

# load tibbles
pvals <- read_tsv( 'pvals.txt.gz', show_col_types = FALSE )
# and true causal info, for AUC
load( 'simtrait.RData' )

# behavior depends on the presence of a true kinship matrix, which tells us if this is a simulation or a real dataset.
# present => sim; absent => real.
is_sim <- file.exists( 'kinship/true.grm.bin' )

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
stopifnot( ncol( pvals ) == n_methods )
aucs <- vector( 'numeric', n_methods )
# names as they appear on the plot
# NOTE: no PCA/LMM marks (will be added to plot separately)
names( aucs ) <- rep.int( kinship_methods$nice, 2 )
# calculate AUCs in desired order
for ( i in 1 : n_methods ) {
    method_code <- method_codes[ i ]
    aucs[i] <- pval_aucpr( pvals[[ method_code ]], causal_indexes)
}

# now make plot of data
dims <- fig_scale( if ( is_sim ) 1.5 else 2 ) # w/h
fig_start(
    'auc',
    width = dims[1],
    height = dims[2],
    mar_l = if ( is_sim ) 11 else 10
)
## boxplot(
##     aucs,
##     horizontal = TRUE,
##     las = 1,
##     xlab = expression(bold(AUC[PR]))
## )
ys <- barplot(
    rev( aucs ),
    horiz = TRUE,
    las = 1,
    xlab = expression(bold(AUC[PR]))
)
# add separating lines
x_line0 <- if ( is_sim ) 0.22 else 0.055
x_line <- -c(x_line0, x_line0)
y_line1 <- ys[ c(1, n_kinship) ]
y_line2 <- ys[ c(1, n_kinship) + n_kinship ]
lines( x_line, y_line1, xpd = NA )
lines( x_line, y_line2, xpd = NA )
# and labels
text( - x_line0 * 1.05, mean( y_line2 ), 'PCA', xpd = NA, srt = 90 )
text( - x_line0 * 1.05, mean( y_line1 ), 'LMM', xpd = NA, srt = 90 )
fig_end()
