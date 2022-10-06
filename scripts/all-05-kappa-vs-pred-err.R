library(readr)
library(ochoalabtools)
library(dplyr)

# a sort of constant
width <- fig_width()

datasets <- c('sim-admix-n1000-m100000-k3-f0.3-s0.5-mc100-h0.8-g20-fes', 'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01')
dataset_names <- c('Admixed Family sim.', '1000 Genomes')

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cccc' )
# we'll only use this to map names for popkin and standard, so exclude WG upfront
kinship_methods <- kinship_methods[ grep( '^wg_', kinship_methods$code, invert = TRUE ), ]
# for convenience, reorder so true is first, then popkins, then standards are last
code_order <- c('true', 'popkin_rom', 'popkin_mor', 'std_rom_lim', 'std_rom', 'std_mor')
kinship_methods <- kinship_methods[ match( code_order, kinship_methods$code ), ]
# now assign colors to each method for plot
kinship_methods$col = 1 : nrow( kinship_methods )

# store processed data from all datasets
data_all <- vector( 'list', length( datasets ) )

# go where the data is
# load precomputed sigmas
setwd( '../data/' )

# navigate datasets
for ( i in 1L : length( datasets ) ) {
    dataset <- datasets[i]
    setwd( dataset )
    data <- read_tsv( 'eigen.txt.gz', col_types = 'cicdd' )
    eval <- read_tsv( 'eval.txt.gz', col_types = 'cid' )
    reml <- read_tsv( 'reml.txt.gz', col_types = 'icdddddddd' )

    # work to correlate V condition numbers to issues with SRMSDs and AUCs
    # to do it broadly (before focusing on particularly troublesome cases), focus on broad cases where this matters
    # apply filters globally (will not use the removed data anymore)
    # first, only V matters according to theory
    data <- data[ data$type == 'V', ]
    # and in evaluations, only LMMs are affected by this argument (and no PCA issues were observed)
    # this is tricky only because the LMM/PCA status is encoded as a substring, not as a separate column
    eval <- eval[ grep( '^lmm_', eval$kinship ), ]
    # and remove that substring now that it's redundant and complicates matching
    eval$kinship <- sub( '^lmm_', '', eval$kinship )

    # sadly they're not aligned, let's align by matching kinship and rep
    data$code <- paste0( data$kinship, ':', data$rep )
    eval$code <- paste0( eval$kinship, ':', eval$rep )
    # reorder!
    eval <- eval[ match( data$code, eval$code ), ]
    # ensure that they're aligned now
    stopifnot( all( data$code == eval$code ) )

    # for clarity remove extra columns now, though not strictly needed
    data$code <- NULL
    data$type <- NULL
    eval$code <- NULL

    # transfer rmsd and auc from eval to data
    data$rmsd <- eval$rmsd
    data$auc <- eval$auc
    
    # time to choose a reference dataset!
    # WG actually has the lowest condition numbers and is more consistent, though STD is a close second
    # separate the three as needed
    data_wg <- data[ grep( '^wg_', data$kinship ), ]
    data_pop <- data[ grep( '^(true|popkin_)', data$kinship ), ]
    data_std <- data[ grep( '^std_', data$kinship ), ]
    # ensure alignment again by creating type columns
    data_wg$type <- sub( '^wg_', '', data_wg$kinship )
    data_std$type <- sub( '^std_', '', data_std$kinship )
    data_pop$type <- sub( '^popkin_', '', data_pop$kinship )
    # handle true correctly
    data_pop$type[ data_pop$type == 'true' ] <- 'rom_lim'
    # and test by matching both type and rep
    data_wg$code <- paste0( data_wg$type, ':', data_wg$rep )
    data_std$code <- paste0( data_std$type, ':', data_std$rep )
    data_pop$code <- paste0( data_pop$type, ':', data_pop$rep )
    # data was already aligned, this confirms it!
    stopifnot( all( data_wg$code == data_std$code ) )
    stopifnot( all( data_wg$code == data_pop$code ) )
    # replace rmsd and auc with deltas from WG (the reference)
    data_pop$auc <- abs( data_pop$auc - data_wg$auc )
    data_pop$rmsd <- abs( data_pop$rmsd - data_wg$rmsd )
    data_std$auc <- abs( data_std$auc - data_wg$auc )
    data_std$rmsd <- abs( data_std$rmsd - data_wg$rmsd )

    # also want the sigma prediction errors, include now!
    # we could include the residual errors but those are much smaller than the genetic sigma errors
    # the marginals suggest there won't be a correlation, but it'll be nice to show that directly
    # calculate the errors as before except using WG as reference, so for popkin that's flipped a bit
    #reml$err_g_std = reml$g_std - reml$g_tr * reml$c_std
    reml$err_g_pop = abs( reml$g_tr - reml$g_wg / reml$c_wg )
    reml$err_g_std = abs( reml$g_std - reml$g_wg * ( reml$c_std / reml$c_wg ) )
    # and add code to match those used earlier
    reml$code <- paste0( reml$type, ':', reml$rep )
    # check alignment now (is already aligned, just confirm!)
    stopifnot( all( reml$code == data_pop$code ) )
    stopifnot( all( reml$code == data_std$code ) )
    # transfer values now!
    data_pop$err_g <- reml$err_g_pop
    data_std$err_g <- reml$err_g_std

    # merge back before replotting?
    data <- bind_rows( data_pop, data_std )
    # use same colors as legend
    data$col <- kinship_methods$col[ match( data$kinship, kinship_methods$code ) ]

    # store processed data
    data_all[[ i ]] <- data
    
    # go back down when done
    setwd( '..' )
}

# start figure vs condition numbers
fig_start(
    'kappa-vs-pred-err',
    width = width,
    height = width * 3/2,
    mar_t = 2
)
# 6 panels in total
par( mfcol = c(3, 2) )

# navigate datasets
for ( i in 1L : length( data_all ) ) {
    data <- data_all[[ i ]]
    
    # plot data of interest!  Wow beautiful correlation!
    plot(
        data$kappa,
        data$err_g,
        col = data$col,
        log = 'x',
        xlab = 'Condition number',
        ylab = 'Sigma_g^2 abs. diff. from WG pred',
        main = dataset_names[i]
    )
    panel_letter( toupper( letters[i] ) ) # , adj = -0.15
    # add legend first time only
    if ( i == 1L )
        legend(
            'topright',
            kinship_methods$nice,
            col = kinship_methods$col,
            pch = 1
        )
    
    plot(
        data$kappa,
        data$rmsd,
        col = data$col,
        log = 'x',
        xlab = 'Condition number',
        ylab = 'RMSD abs. diff. from WG'
    )

    plot(
        data$kappa,
        data$auc,
        col = data$col,
        log = 'x',
        xlab = 'Condition number',
        ylab = 'AUC abs. diff. from WG'
    )
}

fig_end()



# start figure vs sigma error
fig_start(
    'reml-err-vs-pred-err',
    width = width,
    height = width,
    mar_t = 2
)
# 4 panels in total
par( mfcol = c(2, 2) )

# navigate datasets
for ( i in 1L : length( data_all ) ) {
    data <- data_all[[ i ]]
    
    # plot data of interest!  Wow beautiful correlation!
    plot(
        data$err_g,
        data$rmsd,
        col = data$col,
        xlab = 'Sigma_g^2 abs. diff. from WG pred',
        ylab = 'RMSD abs. diff. from WG',
        main = dataset_names[i]
    )
    panel_letter( toupper( letters[i] ) ) # , adj = -0.15
    # add legend first time only
    if ( i == 1L )
        legend(
            'topright',
            kinship_methods$nice,
            col = kinship_methods$col,
            pch = 1
        )
    
    plot(
        data$err_g,
        data$auc,
        col = data$col,
        xlab = 'Sigma_g^2 abs. diff. from WG pred',
        ylab = 'AUC abs. diff. from WG'
    )
}

fig_end()

