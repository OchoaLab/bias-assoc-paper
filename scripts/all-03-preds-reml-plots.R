library(readr)
library(ochoalabtools)
library(popkin)

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cccc' )

# go where the data is
# load precalculated stats
setwd( '../data/' )

# simulation first, use first replicate only
setwd( 'sim-admix-n1000-m100000-k3-f0.3-s0.5-mc100-h0.8-g20-fes' )
data_sim <- read_tsv( 'preds-reml.txt.gz', col_types = 'icdddddddddd' )
setwd( '..' )

# now real data
setwd( 'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01' )
data_real <- read_tsv( 'preds-reml.txt.gz', col_types = 'icdddddddddd' )
setwd( '..' )

boxplot_errors <- function ( data, main = '' ) {
    # types vary per dataset, so determine them here
    types <- unique( data$type )
    # and by component/bias types (constants)
    components <- c( 'e', 'g' )
    short_names <- c( 'st', 'wg' )
    errors <- list()
    labs_type <- c() # grow it in loop
    labs_short <- c() # ditto
    labs_component <- c() # ditto
    # this loop determines boxplot order
    for ( component in components ) {
        for ( type in types ) {
            for ( short_name in short_names ) {
                # name as it appears in column of `data`
                col_name <- paste0( 'err_', component, '_', short_name )
                # separate data as desired
                errors[[ paste0( col_name, '_', type ) ]] <- data[[ col_name ]][ data$type == type ]
                # add pretty names for plot for these properties
                labs_component <- c( labs_component, if ( component == 'e' ) 'Residual' else 'Genetic' )
                labs_type <- c( labs_type, kinship_methods$type[ kinship_methods$code == paste0( 'std_', type ) ] )
                # this is how it matches my other tables
                if ( short_name == 'st' ) short_name <- 'std'
                labs_short <- c( labs_short, kinship_methods$short[ kinship_methods$code == paste0( short_name, '_rom' ) ] )
            }
        }
    }
    # make boxplot!
    boxplot(
        errors,
        main = main,
        names = NA, # individual labels will be plotted with rest
        xaxt = 'n',
        ylab = expression(bold((sigma*minute)^2 - c*sigma^2))
    )
    # reuse popkin-style labeling!  Great for hierarchical/factorial setup
    # though popkin has defaults for these, the raw function doesn't have defaults!
    popkin:::print_labels_multi(
                 labs = cbind( labs_short, labs_type, labs_component ),
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
                 xc_ind = 1 : length( errors ),
                 doMat = FALSE
             )
}

boxplot_emins <- function ( data, main = '' ) {
    # types vary per dataset, so determine them here
    types <- unique( data$type )
    # and by matrix/bias types (constants)
    short_names <- c( 'tr', 'st', 'wg' )
    matrices <- c('k', 'v')
    emins <- list()
    labs_type <- c() # grow it in loop
    labs_short <- c() # ditto
    labs_matrix <- c() # ditto
    # this loop determines boxplot order
    for ( matrice in matrices ) { # fancyful singular to avoid "matrix"
        for ( type in types ) {
            for ( short_name in short_names ) {
                # name as it appears in column of `data`
                col_name <- paste0( 'emin_', matrice, '_', short_name )
                # separate data as desired
                emins[[ paste0( col_name, '_', type ) ]] <- data[[ col_name ]][ data$type == type ]
                # add pretty names for plot for these properties
                labs_matrix <- c( labs_matrix, if ( matrice == 'k' ) 'Kinship' else 'Trait Cov.' )
                labs_type <- c( labs_type, kinship_methods$type[ kinship_methods$code == paste0( 'std_', type ) ] )
                # this is how it matches my other tables
                code <- short_name
                if ( short_name == 'st' ) code <- 'std'
                if ( short_name == 'tr' ) {
                    code <- if ( type == 'rom_lim' ) 'true' else 'popkin_rom'
                } else
                    code <- paste0( code, '_rom' )
                labs_short <- c( labs_short, kinship_methods$short[ kinship_methods$code == code ] )
            }
        }
    }
    # make boxplot!
    boxplot(
        emins,
        main = main,
        names = NA, # individual labels will be plotted with rest
        xaxt = 'n',
        ylab = 'Min. eigenvalue'
    )
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
                 xc_ind = 1 : length( emins ),
                 doMat = FALSE
             )
}


width <- fig_width()
fig_start(
    'preds-reml-errors',
    width = width,
    height = width / 2,
    mar_t = 1.5,
    mar_b = 7
)
par( mfrow = c(1,2) )
par( oma = c(1, 0, 0, 0) )
boxplot_errors( data_sim, main = 'Admixed Family sim.' )
panel_letter('A', adj = -0.2)
boxplot_errors( data_real, main = '1000 Genomes' )
panel_letter('B', adj = -0.2)
mtext( 'Variance Component, Kinship Estimate', side = 1, outer = TRUE )
fig_end()

# interpretation
# - outliers are only real g ROM, and for both bias types the outliers are close, so popkin ROM must be the problem
#   - those outliers are very negative, meaning that popkin ROM's sigma is overestimated (too large) vs theory

fig_start(
    'preds-reml-emins',
    width = width,
    height = width / 2,
    mar_t = 1.5,
    mar_b = 7
)
par( mfrow = c(1,2) )
par( oma = c(1, 0, 0, 0) )
boxplot_emins( data_sim, main = 'Admixed Family sim.' )
panel_letter('A', adj = -0.2)
boxplot_emins( data_real, main = '1000 Genomes' )
panel_letter('B', adj = -0.2)
mtext( 'Matrix type, Kinship Estimate', side = 1, outer = TRUE )
fig_end()

# so the story is consistent that posdef breakdowns for V largely depend on K's status
# why are MOR estimates bad?  for popkin they are not in the guaranteed form!  MOR form exacerbates that, appears ROM is more robust to that issue
