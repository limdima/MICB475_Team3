# Meeting Agenda
1. Proposal grading question: double marks docked for not describing the dataset in detail.
2. Does it matter which condition is set as a reference for DESeq analysis? (DESeq2.R lines 18, 23, 28)
# Meeting Notes
* Grading of the proposal will be discussed with Evelyn.
* We need to figure out how to filter the rooted tree. Hans will go back to the code and look into it. There may be a function that we can use, like “prune samples” that may be used for rooted trees. The prune function may also allow us to filter after making the phyloseq object. On this, please create a directory in GitHub to add our code.
  *  One option is to filter our data first before creating any of the objects like the tree. We have a couple options: (1) Pruning our data with qiime2, (2) Filtering in R and reimporting into qiime2 to create the various phyloseq components. As long as get 515 samples by the end.
* From a Venn Diagram created for the Core Microbiome analysis, there were no unique genera to only PPMS or SPMS. RPMS had some unique taxa, and most were shared between all three. This difference between PPMS/SPMS and RPMS might be of interest to us. It would also be interesting to compare this to patients with no MS.
  * Either we add another circle to the Venn diagram, or remove non-unique genera between any MS patients and controls. We should also increase our prevalence treshhold, even if that gets rid of some of our shared genera between PPMS and SPMS.
  * We should investigate the four genera that are unique to SPMS and PPMS.
* It would be good to compare with the control to analyse differential abundance. ANCOM or ALDEX2 (both plugins for qiime2) might be better for our purposes, since DESEQ is not made by 16S data.
* Aim 1 gave expected data, no significant differences between treated and untreated, and significant beta diversity but on a very small scale.
