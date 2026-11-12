Baseline scripts (reference copies, unchanged)
==============================================

These are the original scripts the StRym analysis steps were generalized from. They were written for one experiment
(soybean heat response: genotypes x time points, "heat" vs "control") and read a folder of per-genotype, per-time-point
Cuffdiff files. StRym does NOT run these files; it runs the generalized versions in the parent folder:

  mse_simple.R           ->  ../modified_shannon_entropy.R   (+ shared code in ../common.R)
  shannon_simple.R       ->  ../shannon_entropy.R            (+ shared code in ../common.R)
  pca_plot_script.R      ->  ../pca.R
  simple_trajectories.R  ->  ../trajectory.R
  single_cell/*.Rmd      ->  ../seurat_clustering.R
  DESeq2 (link only: github.com/zhangdongqin/easy_deseq2, script deseq.r)  ->  ../deseq2.R  (explicit contrasts; the original's
                         results() call has no contrast and it runs human GO / KEGG enrichment, which is left out)

They are kept so the generalized code can be checked against them: tests/testthat/test-baselines.R runs the original
mse_simple.R and shannon_simple.R on a synthetic experiment and requires StRym's results to be identical.
