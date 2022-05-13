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
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double"),
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep <- opt$rep
herit <- opt$herit
m_causal_fac <- opt$m_causal_fac
fes <- opt$fes

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
setwd( '../data/' )
setwd( name )

# check that replicate doesn't already exist (don't overwrite!)
dir_rep <- paste0( 'rep-', rep )
if ( dir.exists( dir_rep ) )
    stop( 'Replicate ', rep, ' already exists!' )
dir.create( dir_rep )

################
### simtrait ###
################

# load precalculated popkin kinship matrix
kinship <- read_grm( 'kinship/popkin_rom' )$kinship

# load FAM table, for phen output
fam <- read_fam( name_in )

# now that we have number of individuals, use a tenth of that for m_causal, rounding
m_causal <- round( nrow( fam ) / m_causal_fac )
message( 'm_causal: ', m_causal )

# load genotype data with BEDMatrix
X <- BEDMatrix( name_in, simple_names = TRUE )

# simulate trait
obj_trait <- sim_trait(
    X = X,
    m_causal = m_causal,
    herit = herit,
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

# move into replicate now, where output data is stored
setwd( dir_rep )

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
