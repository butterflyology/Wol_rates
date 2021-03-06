### Repo for Marion and Hamm (2016) "A hierarchical Bayesian approach to estimate endosymbiont infection rates" for a special issue of in *Frontiers in Microbiology*
Release 1.0 of code. [![DOI](https://zenodo.org/badge/51758300.svg)](https://zenodo.org/badge/latestdoi/51758300)

Major commits by date (2016):

  * 15 February - First commit
  * 11 April - First calculatoin of vcv matrix based on phylogeny
  * 20 April - Corrected Erebidae naming error in Regier et al. (2013) data
  * 24 April - Added OU correction with multiple alpha value
  * 25 April - Phylogenetically corrected model, early code for simulating data
  * 26 April - Models incorporating phylogenetic correction working, full simuation code working for BM and OU models
  * 27 April - Made tree ultrametric, calculated vcv from ultrametric tree, ran models with different vcv corrections
  * 28 April - Made figures, barplots by samples and proportion of species in family, added updated .tex files.
  * 29 April - Updated stan files and modified plots
  * 4 May - Updated plots, added model averaging
  * 8 May - Updated ms latex files
  * 14 May - Corrected model averaging error for family level predictions
  * 24 May - Added correlation and plots for association between CI and sample size
  * 14 June - Updated repo to reflect state upon submission of ms. Latex docs in folder, code in submitted state. Fingers crossed. 
  * 20 June - Removed the paper files, too soon to have them posted on git (should have been on ArXiv, honestly)
  * 11 October - Ammended code to accomodate changes in stan, added code for posterior predictive checks of the model.
  * 1 November - Ammended code to work with new "yaRrr" package
  * 15 November - Release version of cade created.

Major commits by date (2017):
  * 21 January - Updated code to use `spaceMovie` color palette, started issue to correct depreciate (but still functional) code for `stan`.
