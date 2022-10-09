library(readr)
library(ochoalabtools)
library(popkin) # for :::print_labels_multi

# threshold for negative eigenvalues
cut_evs <- -1e-7
# a sort of constant
width <- fig_width()

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cccc' )

# go where the data is
# load precomputed sigmas
setwd( '../data/' )

# simulation first, use first replicate only
setwd( 'sim-admix-n1000-m100000-k3-f0.3-s0.5-mc100-h0.8-g20-fes' )
data_sim <- read_tsv( 'eigen.txt.gz', col_types = 'cicdid' )
eval_sim <- read_tsv( 'eval.txt.gz', col_types = 'cid' )
setwd( '..' )

# now real data
setwd( 'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01' )
data_real <- read_tsv( 'eigen.txt.gz', col_types = 'cicdid' )
eval_real <- read_tsv( 'eval.txt.gz', col_types = 'cid' )
setwd( '..' )


boxplot_data <- function ( data, main = '', is_sim = TRUE, col_name = 'emin', ylab = 'Min. eigenvalue', cut = NA, ... ) {
    # for version of analysis that depends on threshold, apply threshold now, overwrite column (locally only)
    if ( !is.na( cut ) )
        data[[ col_name ]] <- data[[ col_name ]] < cut
    
    types <- c('rom_lim', 'rom', 'mor')
    # no true/limits for real data
    if ( !is_sim )
        types <- types[-1]
    short_names <- c( 'popkin', 'std', 'wg' )
    matrices <- c('kinship', 'V')
    data2 <- list()
    labs_type <- c() # grow it in loop
    labs_short <- c() # ditto
    labs_matrix <- c() # ditto
    # this loop determines boxplot order
    for ( matrice in matrices ) { # fancyful singular to avoid "matrix"
        for ( type in types ) {
            for ( short_name in short_names ) {
                # get name that matches "kinship" column in input data
                kinship_name <- paste0( short_name, '_', type )
                # map only weird case that doesn't fit pattern
                if ( is_sim && kinship_name == 'popkin_rom_lim' )
                    kinship_name <- 'true'
                # create unique output names, though exactly what they are doesn't matter
                out_name <- paste0( kinship_name, '_', matrice )
                # separate data as desired
                data2[[ out_name ]] <- data[[ col_name ]][ data$kinship == kinship_name & data$type == matrice ]
                # add pretty names for plot for these properties
                labs_matrix <- c( labs_matrix, if ( matrice == 'kinship' ) 'Kinship' else 'Trait Cov.' )
                labs_type <- c( labs_type, kinship_methods$type[ kinship_methods$code == paste0( 'std_', type ) ] )
                labs_short <- c( labs_short, kinship_methods$short[ kinship_methods$code == kinship_name ] )
            }
        }
    }
    if ( is.na( cut ) ) {
        # make boxplot!
        boxplot(
            data2,
            main = main,
            names = NA, # individual labels will be plotted with rest
            xaxt = 'n',
            ylab = ylab,
            ...
        )
        # location of boxes
        xc_ind <- 1 : length( data2 )
    } else {
        # thresholded analysis is better as barplot
        # bar locations here are non-trivial
        xc_ind <- barplot(
            sapply( data2, mean ),
            main = main,
            names = NA, # individual labels will be plotted with rest
            xaxt = 'n',
            ylab = ylab,
            ...
        )
    }
    # reuse popkin-style labeling!  Great for hierarchical/factorial setup
    # though popkin has defaults for these, the raw function doesn't have defaults!
    popkin:::print_labels_multi(
                 labs = cbind( labs_short, labs_type, labs_matrix ),
                 labs_cex = c(0.7, 0.6, 1),
                 labs_las = c(2, 0, 0),
                 labs_line = c(0.5, 4, 6),
                 labs_lwd = 1, # default
                 labs_sep = c(FALSE, TRUE, TRUE),
                 labs_even = c(FALSE, TRUE, TRUE),
                 labs_ticks = FALSE, # default
                 labs_text = TRUE, # default
                 labs_col = 'black', # default
                 # these align with barplot-specific changes
                 xc_ind = xc_ind,
                 doMat = FALSE
             )
}

# plot minimum eigenvalues!
fig_start(
    'emin',
    width = width,
    height = width / 2,
    mar_t = 1.5,
    mar_b = 7
)
par( mfrow = c(1,2) )
par( oma = c(1, 0, 0, 0) )
boxplot_data( data_sim, main = 'Admixed Family sim.' )
panel_letter('A', adj = -0.2)
boxplot_data( data_real, main = '1000 Genomes', is_sim = FALSE )
panel_letter('B', adj = -0.2)
mtext( 'Matrix type, Kinship Estimate', side = 1, outer = TRUE )
fig_end()

# version that applies threshold to report non-posdef proportions
fig_start(
    'emin-cut',
    width = width,
    height = width / 2,
    mar_t = 1.5,
    mar_b = 7
)
par( mfrow = c(1,2) )
par( oma = c(1, 0, 0, 0) )
boxplot_data( data_sim, main = 'Admixed Family sim.', cut = cut_evs, ylab = 'Proportion non-posdef' )
panel_letter('A', adj = -0.2)
boxplot_data( data_real, main = '1000 Genomes', is_sim = FALSE, cut = cut_evs, ylab = 'Proportion non-posdef' )
panel_letter('B', adj = -0.2)
mtext( 'Matrix type, Kinship Estimate', side = 1, outer = TRUE )
fig_end()

# plot condition numbers
# this plot looks way better with a shared y axis across panels
ylim = range( data_sim$kappa, data_real$kappa )
fig_start(
    'kappa',
    width = width,
    height = width / 2,
    mar_t = 1.5,
    mar_b = 7
)
par( mfrow = c(1,2) )
par( oma = c(1, 0, 0, 0) )
boxplot_data( data_sim, main = 'Admixed Family sim.', col_name = 'kappa', ylab = 'Condition number', log = 'y', cex.axis = 0.7, ylim = ylim )
panel_letter('A', adj = -0.2)
boxplot_data( data_real, main = '1000 Genomes', is_sim = FALSE, col_name = 'kappa', ylab = 'Condition number', log = 'y', cex.axis = 0.7, ylim = ylim )
panel_letter('B', adj = -0.2)
mtext( 'Matrix type, Kinship Estimate', side = 1, outer = TRUE )
fig_end()

# just add test for number of negative eigenvectors never exceeding 1...
stopifnot( max( data_sim$n_neg_evs ) <= 1 )
stopifnot( max( data_real$n_neg_evs ) <= 1 )

# so the story is consistent that posdef breakdowns for V largely depend on K's status
# why are MOR estimates bad?  for popkin they are not in the guaranteed form!  MOR form exacerbates that, appears ROM is more robust to that issue

