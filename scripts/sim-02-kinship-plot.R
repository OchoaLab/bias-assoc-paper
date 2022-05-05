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
setwd( 'kinship' )

# kinship files have these names, will appear in this order
codes <- c(
    'true',
    'popkin_rom',
    'popkin_mor',
    'std_rom_lim',
    'std_rom',
    'std_mor',
    'gcta_rom_lim',
    'gcta_mor',
    'wg_rom_lim',
    'wg_rom',
    'wg_mor'
)
# map to human-readable names from table
titles <- kinship_methods$nice[ match( codes, kinship_methods$code ) ]

# read in all kinship matrices
data <- lapply( codes, function ( name ) {
    # all are 2x, scaled as GCTA wants them, halve here for plots
    read_grm( name )$kinship / 2
})

# HACK, since we don't have gcta_rom yet, to have a nice figure with columns aligned, insert a NULL in the right place
data <- c( data[ 1:7 ], list(NULL), data[ 8:length(data) ] )

# save figure in lower level
setwd( '..' )

# visualize all matrices for test
dims <- fig_scale( ratio = 3/4 )
fig_start(
    'kinship',
    width = dims[1],
    height = dims[2]
)
plot_popkin(
    inbr_diag( data ),
    titles = titles,
    layout_rows = 4,
    mar = c(0, 2),
    panel_letters_adj = 0 # old default, works better here because there's no labels
)
fig_end()
