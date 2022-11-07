# install CRAN packages (include only packages not already installed)
install.packages( c('devtools', 'tibble', 'readr', 'dplyr', 'optparse', 'BEDMatrix', 'genio', 'popkin', 'bnpsd', 'simfam', 'simtrait', 'Matrix') )

# github-only packages
library(devtools)
install_github( c('OchoaLab/popkinsuppl', 'OchoaLab/ochoalabtools', 'OchoaLab/genbin', 'OchoaLab/simgenphen') )
