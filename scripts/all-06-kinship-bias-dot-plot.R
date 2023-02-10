library(genio) # read_grm
library(popkin) # inbr
library(ochoalabtools) # fig_start

# go where the data is
setwd( '../data/' )

# start figure in this location
scale <- fig_scale( 2 ) # w/h ratio
fig_start(
    'kinship-bias',
    width = scale[1],
    height = scale[2],
    mar_t = 1.5
)
# two panels
par( mfrow = c(1, 2) )

# simulation first
# just use the same old simulation, primary of paper
setwd( 'sim-admix-n1000-m100000-k3-f0.3-s0.5-mc100-h0.8-g20-fes' )

# use first replicate only
setwd( 'rep-1/kinship' )

# load the three bias types kinship matrices, limits only
# transform diagonals to inbreeding immediately
kinship_true <- inbr_diag( read_grm( 'true' )$kinship / 2 )
kinship_std_rom_lim <- inbr_diag( read_grm( 'std_rom_lim' )$kinship / 2 )
kinship_wg_rom_lim <- inbr_diag( read_grm( 'wg_rom_lim' )$kinship / 2 )
# set range of plot, symmetric so its more pleasing
kin_range <- range( kinship_true, kinship_std_rom_lim, kinship_wg_rom_lim )

# initialize plot area
plot(
    NA,
    xlim = kin_range,
    ylim = kin_range,
    xlab = 'True kinship',
    ylab = 'Biased ROM lim.',
    main = 'Admixed Family sim.'
)
panel_letter( 'A', adj = -0.2 )
# plot y=x diagonal in background
abline( 0, 1, lty = 2, col = 'grey' )

# plot off-diagonal without duplicating dots
indexes <- lower.tri( kinship_true )
points( kinship_true[ indexes ], kinship_std_rom_lim[ indexes ], pch = '.', col = 'red' )
# diagonal use diff color
points( diag( kinship_true ), diag( kinship_std_rom_lim ), pch = '.', col = 'orange' )

# repeat WG (will be on top of STD)
points( kinship_true[ indexes ], kinship_wg_rom_lim[ indexes ], pch = '.' )
points( diag( kinship_true ), diag( kinship_wg_rom_lim ), pch = '.', col = 'blue' )

legend(
    'topleft',
    c( 'Standard ROM lim. kinship', 'Standard ROM lim. inbreeding', 'WG ROM lim. kinship', 'WG ROM lim. inbreeding' ),
    col = c( 'red', 'orange', 'black', 'blue' ),
    pch = 19,
    bty = 'n',
    cex = 0.7
)


# now real data
setwd( '../../../tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01/kinship' )

# load the three bias types kinship matrix estimates (no "true" for real data), ROM estimates only so popkin's is the unbiased one
kinship_popkin_rom <- inbr_diag( read_grm( 'popkin_rom' )$kinship / 2 )
kinship_std_rom <- inbr_diag( read_grm( 'std_rom' )$kinship / 2 )
kinship_wg_rom <- inbr_diag( read_grm( 'wg_rom' )$kinship / 2 )

# set range of plot, symmetric so its more pleasing
kin_range <- range( kinship_popkin_rom, kinship_std_rom, kinship_wg_rom )

# initialize plot area
plot(
    NA,
    xlim = kin_range,
    ylim = kin_range,
    xlab = 'Popkin ROM est.',
    ylab = 'Biased ROM est.',
    main = '1000 Genomes'
)
panel_letter( 'B', adj = -0.2 )
# plot y=x diagonal in background
abline( 0, 1, lty = 2, col = 'grey' )

# plot off-diagonal without duplicating dots
indexes <- lower.tri( kinship_popkin_rom )
points( kinship_popkin_rom[ indexes ], kinship_std_rom[ indexes ], pch = '.', col = 'red' )
# diagonal use diff color?
points( diag( kinship_popkin_rom ), diag( kinship_std_rom ), pch = '.', col = 'orange' )

# repeat WG (will be on top of STD)
points( kinship_popkin_rom[ indexes ], kinship_wg_rom[ indexes ], pch = '.' )
points( diag( kinship_popkin_rom ), diag( kinship_wg_rom ), pch = '.', col = 'blue' )

legend(
    'topleft',
    c( 'Standard ROM est. kinship', 'Standard ROM est. inbreeding', 'WG ROM est. kinship', 'WG ROM est. inbreeding' ),
    col = c( 'red', 'orange', 'black', 'blue' ),
    pch = 19,
    bty = 'n',
    cex = 0.7
)

fig_end()

