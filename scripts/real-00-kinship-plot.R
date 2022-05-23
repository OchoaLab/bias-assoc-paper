library(optparse) # for terminal options
library(genio)    # to read GRM files
library(popkin)   # to plot
library(ochoalabtools) # for nice PDF
library(readr)

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

# real data has the issue or sample ordering!
# read FAM for current order and subpop correspondence
fam <- read_fam( 'data' )
# read annotations
subpop_info <- read_tsv( 'pops-annot.txt', comment = '#', show_col_types = FALSE )
# map subpopulations using sub-subpopulations
fam$superpop <- subpop_info$superpop[ match( fam$fam, subpop_info$pop ) ]
# reorder individuals so subpopulations come out in desired order:
indexes <- order( match( fam$fam, subpop_info$pop ) )
# apply reordering to fam (will do kinship matrices as they are loaded)
fam <- fam[ indexes, ]

setwd( 'kinship' )

# kinship files have these names, will appear in this order
codes <- c(
    'popkin_rom',
    'popkin_mor',
    'std_rom',
    'std_mor',
    'wg_rom',
    'wg_mor'
)
# map to human-readable names from table
titles <- kinship_methods$nice[ match( codes, kinship_methods$code ) ]

# read in all kinship matrices
data <- lapply( codes, function ( name ) {
    # all are 2x, scaled as GCTA wants them, halve here for plots
    kinship <- read_grm( name )$kinship / 2
    # apply reordering to everything!
    kinship <- kinship[ indexes, indexes ]
    # also transform diagonal, return
    return( inbr_diag( kinship ) )
})

# the popkin MOR kinship matrix generally has higher values, we want to preserve that range
# the biased kinship matrice have usually lower values but the diagonal has much larger values, those are the ones we want to cap
alpha <- 0.01
kinship_max <- quantile( diag( data[[2]] ), probs = 1 - alpha )
#kinship_max <- max( data[[1]] )
cap_kinship <- function( x ) {
    x[ x > kinship_max ] <- kinship_max
    return( x )
}
data <- lapply( data, cap_kinship )

# save figure in lower level
setwd( '..' )

# visualize all matrices for test
dims <- fig_scale( ratio = 2.2/3 ) # w/h
fig_start(
    'kinship',
    width = dims[1],
    height = dims[2] # hack
)
plot_popkin(
    data,
    titles = titles,
    labs = fam$superpop,
    labs_line = 0.2,
    labs_cex = 0.7,
    layout_rows = 3,
    mar = c(2, 2)
#    panel_letters_adj = 0 # old default, works better here because there's no labels
)
fig_end()
