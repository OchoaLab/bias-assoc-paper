# this scripts simulates a random trait for the given real dataset, storing key values

library(optparse)
library(simtrait)
library(genio)
library(BEDMatrix)

# the name is for dir only, actual file is just "data"
name_in <- 'data'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option(c("-r", "--rep"), type = "integer", default = 1, 
                help = "Replicate number", metavar = "int"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep <- opt$rep
herit <- opt$herit
fes <- opt$fes

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
setwd( '../data/' )
setwd( name )

# if this directory exists here, this is real data
is_sim <- !file.exists( 'kinship' )

# check that replicate doesn't already exist (don't overwrite!)
dir_rep <- paste0( 'rep-', rep )
# create directory if needed
if ( !dir.exists( dir_rep ) )
    dir.create( dir_rep )

# if it's a simulation, enter replicate because genotypes and kinship matrices are in there
# otherwise stay in current directory
if ( is_sim )
    setwd( dir_rep )

################
### simtrait ###
################

# set missing parameters to null
p_anc <- NULL
kinship <- NULL
# load non-null case as needed
if ( is_sim ) {
    load( 'p_anc.RData' ) # loads p_anc
} else {
    # load precalculated popkin kinship matrix
    kinship <- read_grm( 'kinship/popkin_rom' )$kinship
}

# load FAM table, for phen output
fam <- read_fam( name_in )

# now that we have number of individuals, use it to calculate m_causal, rounding
# always set m_causal automatically using the other variables, as done for pca-assoc project
# for default n,h values it equals the previous default m_causal too!
m_causal <- round( nrow( fam ) * herit / 8 )
message( 'm_causal: ', m_causal )

# load genotype data with BEDMatrix
X <- BEDMatrix( name_in, simple_names = TRUE )

# simulate trait
obj_trait <- sim_trait(
    X = X,
    m_causal = m_causal,
    herit = herit,
    p_anc = p_anc,
    kinship = kinship,
    fes = fes
)
# extract data of interest
# trait vector
trait = obj_trait$trait
# randomly-picked causal locus index
causal_indexes = obj_trait$causal_indexes
# locus coefficient vector (for causal loci only)
causal_coeffs = obj_trait$causal_coeffs

# for real data, move into replicate now, where output data is stored
# (for sims we're alread where we need to be)
if ( !is_sim )
    setwd( dir_rep )

# create and move into new directory if heritability is non-default
if ( herit != 0.8 ) {
    dir_out <- paste0( 'h-', herit )
    if ( !dir.exists( dir_out ) )
        dir.create( dir_out )
    setwd( dir_out )
}

# now save, as R data
save(
    trait,
    causal_indexes,
    causal_coeffs,
    file = 'simtrait.RData'
)

# save trait as phen file too, for GCTA
fam$pheno <- trait
write_phen( name_in, fam )
