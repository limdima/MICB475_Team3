# Notebook

# October 16, 2024

## Removing mitochondria and chloroplast sequences
qiime taxa filter-table \
  --i-table ms-table.qza \
  --i-taxonomy ms-taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table ms-table-no-mitochondria-no-chloroplast.qza

## Visualize mitochondria and chloroplast filtered ASVs stats
qiime feature-table summarize \
  --i-table ms-table-no-mitochondria-no-chloroplast.qza \
  --o-visualization ms-table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/MS/corrected_ms_metadata.tsv

## Frequency-based filtering (ASV's with <0.005% of total reads are filtered out as they may be sequencing errors rather than true biological variants)
qiime feature-table filter-features \
--i-table ms-table-no-mitochondria-no-chloroplast.qza \
--p-min-frequency ___ \
--o-filtered-table ms-mit-chlor-freq-filtered-table.qza


## Visualize mitochondria, chloroplast, and frequency-based filtered ASVs stats
qiime feature-table summarize \
  --i-table ms-mit-chlor-freq-filtered-table.qza \
  --o-visualization ms-mit-chlor-freq-filtered-table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/MS/corrected_ms_metadata.tsv

## Transfering the mitochondria, chloroplast, and frequency-based filtered table to local computer.
scp root@10.19.139.120:/data/project2/ms-mit-chlor-freq-filtered-table.qzv .

## Generating a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ms-rep-seqs.qza \
  --o-alignment ms-aligned-rep-seqs.qza \
  --o-masked-alignment ms-masked-aligned-rep-seqs.qza \
  --o-tree ms-unrooted-tree.qza \
  --o-rooted-tree ms-rooted-tree.qza 

## Alpha-rarefaction
qiime diversity alpha-rarefaction \
  --i-table ms-mit-chlor-freq-filtered-table.qza \
  --i-phylogeny ms-rooted-tree.qza \
  --p-max-depth ____ \
  --m-metadata-file /mnt/datasets/project_2/MS/corrected_ms_metadata.tsv \
  --o-visualization ms-alpha-rarefaction.qzv

# October 15: Dataset Processing for Team Proposal

## Importing and demultiplexing MS data to the project2 directory
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /mnt/datasets/project_2/MS/ms_manifest.tsv \
  --output-path ./ms_demux_seqs.qza

## Create visualization of demultiplexed samples
qiime demux summarize \
  --i-data ms_demux_seqs.qza \
  --o-visualization ms_demux_seqs.qzv

## Copy file to personal computer for visualization
scp root@10.19.139.120:/data/project2/ms_demux_seqs.qzv .

## Determine ASVs with DADA2
#medians at each position are all above 32 so keeping whole read (151bp)
qiime dada2 denoise-single \
  --i-demultiplexed-seqs ms_demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 151 \
  --o-representative-sequences ms-rep-seqs.qza \
  --o-table ms-table.qza \
  --o-denoising-stats ms-stats.qza

##Visualize DADA2 stats
qiime metadata tabulate \
    --m-input-file ms-stats.qza \
    --o-visualization ms-stats.qzv

## Visualize ASVs stats
qiime feature-table summarize \
  --i-table ms-table.qza \
  --o-visualization ms-table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/MS/corrected_ms_metadata.tsv
  
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization ms-rep-seqs.qzv


## Extracting amplicon of interest from the reference database (primer sequences found in Konnert et al. UJEMI paper)
qiime feature-classifier extract-reads \
  --i-sequences /mnt/datasets/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 151 \
  --o-reads ms-ref-seqs-trimmed.qza

## Training classifier with your new ref-seq file
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ms-ref-seqs-trimmed.qza \
  --i-reference-taxonomy /mnt/datasets/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier ms-classifier.qza

## Use the trained classifier to assign taxonomy to our reads
qiime feature-classifier classify-sklearn \
  --i-classifier ms-classifier.qza \
  --i-reads ms-rep-seqs.qza \
  --o-classification ms-taxonomy.qza

## Next Steps: Taxonomic analysis
