
# generate data on standard K=3 admixture
# time Rscript sim-00-mk-kinship-stats.R 
# # 0m28.729s viiiaR5
# # 0m48.462s ideapad
# time Rscript sim-00-mk-kinship-stats.R --fes
# inv trait and family!  This is the only version shown in paper so far!!!
# also only version redone with most recent changes (simfam and genbin)
time Rscript sim-00-mk-kinship-stats.R -g 20 --fes
# 0m53.873s ideapad

# kinship plots
# should be about the same for both trait types, so just run one
# (made both to just load the most convenient copy while we decide)
# time Rscript sim-01-kinship-plot.R
# time Rscript sim-01-kinship-plot.R --fes
time Rscript sim-01-kinship-plot.R -g 20 --fes

# statistic correlation heatmaps
# Rscript sim-02-pval-beta-corr.R
# Rscript sim-02-pval-beta-corr.R --fes
Rscript sim-02-pval-beta-corr.R -g 20 --fes
# Subset: PCA, Weir-Goudet lim., PCA, Standard ROM lim., PCA, GCTA lim., PCA, Weir-Goudet est., PCA, Standard ROM est., PCA, Standard MOR est., PCA, GCTA est.
# Range: 0.995473892620475, 1
# Subset: PCA, Weir-Goudet lim., PCA, Standard ROM lim., PCA, GCTA lim., PCA, Popkin est., PCA, Weir-Goudet est., PCA, Standard ROM est., PCA, Standard MOR est., PCA, GCTA est.
# Range: 0.97312763475995, 1
# Subset: PCA, True Kinship, PCA, Weir-Goudet lim., PCA, Standard ROM lim., PCA, GCTA lim., PCA, Popkin est., PCA, Weir-Goudet est., PCA, Standard ROM est., PCA, Standard MOR est., PCA, GCTA est.
# Range: 0.970466941724987, 1
# Subset: LMM, True Kinship, LMM, Weir-Goudet lim., LMM, Standard ROM lim.
# Range: 0.999999839327728, 1
# Subset: LMM, True Kinship, LMM, Weir-Goudet lim., LMM, Standard ROM lim., LMM, GCTA lim.
# Range: 0.991417155638545, 1
# Subset: LMM, Popkin est., LMM, Weir-Goudet est., LMM, Standard ROM est.
# Range: 0.999997467390203, 1
# Subset: LMM, Standard MOR est., LMM, GCTA est.
# Range: 1, 1
# Subset: LMM, Popkin est., LMM, Weir-Goudet est., LMM, Standard ROM est., LMM, Standard MOR est., LMM, GCTA est.
# Range: 0.967053222476237, 1
# Subset: LMM, True Kinship, LMM, Weir-Goudet lim., LMM, Standard ROM lim., LMM, GCTA lim., LMM, Popkin est., LMM, Weir-Goudet est., LMM, Standard ROM est., LMM, Standard MOR est., LMM, GCTA est.
# Range: 0.877022643496176, 1

# Rscript sim-03-auc.R
# Rscript sim-03-auc.R --fes
Rscript sim-03-auc.R -g 20 --fes


### THEORY/OBSOLETE ###

# empirically test the hypothesis that centering and EVD approximation commute
# this script simulates kinship as usual, we use r = K
# RESULTS:
# - commutability is not exact, especially noticeable under family structure, but still appears to be a very good approximation
# - the intercept is in the rowspace in both cases, as expected
# - rowspaces are a mess (centering does alter rowspaces enough that they really don't match; see last line of each run below, calculated via matrix ranks on concatenations of matrices)
# 
# 1) version without family structure (truly low-dimensional, K=3, Fst=0.3)
Rscript test-00-rank-pca-center.R
# RMSD center/dim-r, vs dim-r/center: 0.000501087785158289
# Corr center/dim-r, vs dim-r/center: 0.999993595141344
# Sum of kinship, center/dim-r: 1.3e-14
# Sum of kinship, dim-r/center: -3.7e-11
# Rank center/dim-r: 3
# Rank dim-r/center: 3
# Rank center/dim-r + dim-r/center: 6
# Rank center/dim-r + dim-r: 6
# Rank dim-r: 3
# Rank dim-r/center: 3
# Rank dim-r/center + dim-r: 4
# Rank intercept + dim-r: 4
# Rank intercept + dim-r/center: 4
# Rank intercept + center/dim-r: 4
# Rank intercept + dim-r/center + dim-r: 4
# Rank intercept + center/dim-r + dim-r: 7
#
# 2) larger K, lower FST, same result
Rscript test-00-rank-pca-center.R -k 10 -f 0.1
# RMSD center/dim-r, vs dim-r/center: 0.0005038297388324
# Corr center/dim-r, vs dim-r/center: 0.999923201852988
# Sum of kinship, center/dim-r: -1.95e-15
# Sum of kinship, dim-r/center: 1.9e-11
# Rank center/dim-r: 10
# Rank dim-r/center: 10
# Rank center/dim-r + dim-r/center: 19
# Rank center/dim-r + dim-r: 20
# Rank dim-r: 10
# Rank dim-r/center: 10
# Rank dim-r/center + dim-r: 11
# Rank intercept + dim-r: 11
# Rank intercept + dim-r/center: 11
# Rank intercept + center/dim-r: 11
# Rank intercept + dim-r/center + dim-r: 11
# Rank intercept + center/dim-r + dim-r: 20
#
# 3) back to default (K=3, Fst=0.3), but add family structure
Rscript test-00-rank-pca-center.R -g 20
# RMSD center/dim-r, vs dim-r/center: 0.00309263509709354
# Corr center/dim-r, vs dim-r/center: 0.99962525703413
# Sum of kinship, center/dim-r: -3.28e-14
# Sum of kinship, dim-r/center: 4.45e-11
# Rank center/dim-r: 3
# Rank dim-r/center: 3
# Rank center/dim-r + dim-r/center: 6
# Rank center/dim-r + dim-r: 6
# Rank dim-r: 3
# Rank dim-r/center: 3
# Rank intercept + dim-r: 4
# Rank intercept + dim-r/center: 4
# Rank intercept + center/dim-r: 4
# Rank intercept + dim-r/center + dim-r: 4
# Rank intercept + center/dim-r + dim-r: 7

