# run LMM assoc in R, with detailed reporting to understand what is happening in reality for non-posdef kinship
# will focus on true vs lims for simplicity

library(optparse)
library(readr)
library(genio)
library(tibble)
library(BEDMatrix)
library(dplyr)

# a name for temporary BED/etc data, under project dir
name <- 'data'
# things to process
codes <- c('true', 'wg_rom_lim', 'std_rom_lim')

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

# go where we want data outputs to be
setwd( '../data' )
setwd( dir_out )

# if this directory exists here, this is real data
is_sim <- !file.exists( 'kinship' )

dir_rep <- paste0( 'rep-', rep )
setwd( dir_rep )

#############
### ASSOC ###
#############

# load shared things we need

# load phenotype
y <- read_phen( name )$pheno

Vs_inv <- lapply( codes, function (code) {
    # load precalculated V matrix
    V <- read_grm( paste0( 'V/', code ) )$kinship
    # invert V, store that only
    return( solve( V ) )
})
names( Vs_inv ) <- codes

# get genotypes loaded sort of
# if real, they'll be in a lower level...
name_geno <- if ( is_sim ) name else paste0( '../', name )
# load with low mem, hoping not to actually do whole genome, just a few examples
X <- BEDMatrix( name_geno, simple_names = TRUE )

# output data to grow
data <- NULL

# for the rest we need covariates!
# scan all data!
m_loci <- ncol( X )
for ( i in 1L : m_loci ) {
    # use BEDMatrix orientation (loci on columns!)
    xi <- X[ , i ]
    
    # set up our minimal covariate/design matrix
    # covariates are columns in this notation!
    Z <- cbind( 1, xi )

    # get coefficients for each code

    for ( code in codes ) {
        # compute a shared factor
        ZV_inv <- crossprod( Z, Vs_inv[[ code ]] )
        # now calculate final coefficients
        beta <- solve( ZV_inv %*% Z ) %*% ( ZV_inv %*% y )
        # add a row to output, nice tody format
        data_in <- tibble(
            index = i,
            kinship = code,
            intercept = beta[1],
            beta = beta[2]
        )
        data <- bind_rows( data, data_in )
    }
}

# save data inside rep
write_tsv( data, 'lmm-intercept-test.txt.gz' )
