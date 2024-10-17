# Notebook

## October 15: Dataset Processing for Team Proposal

#(done) importing and demultiplexing MS data to the project2 directory

qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /mnt/datasets/project_2/MS/ms_manifest.tsv \
  --output-path ./ms_demux_seqs.qza

#(done) Create visualization of demultiplexed samples

qiime demux summarize \
  --i-data ms_demux_seqs.qza \
  --o-visualization ms_demux_seqs.qzv

#(done) Copy file to personal computer for visualization
scp root@10.19.139.120:/data/project2/ms_demux_seqs.qzv .

#(done) Determine ASVs with DADA2
#medians at each position are all above 32 so keeping whole read (151bp)
#start 2:13PM; end 

qiime dada2 denoise-single \
  --i-demultiplexed-seqs ms_demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 151 \
  --o-representative-sequences ms-rep-seqs.qza \
  --o-table ms-table.qza \
  --o-denoising-stats ms-stats.qza

#(done)Visualize DADA2 stats

qiime metadata tabulate \
    --m-input-file ms-stats.qza \
    --o-visualization ms-stats.qzv

#(done) Visualize ASVs stats

qiime feature-table summarize \
  --i-table ms-table.qza \
  --o-visualization ms-table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/MS/corrected_ms_metadata.tsv
  

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization ms-rep-seqs.qzv


#(done) Extracting amplicon of interest from the reference database (primer sequences found in Konnert et al. UJEMI paper)

qiime feature-classifier extract-reads \
  --i-sequences /mnt/datasets/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 151 \
  --o-reads ms-ref-seqs-trimmed.qza

#(done)Training classifier with your new ref-seq file

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ms-ref-seqs-trimmed.qza \
  --i-reference-taxonomy /mnt/datasets/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier ms-classifier.qza

#Use the trained classifier to assign taxonomy to our reads

qiime feature-classifier classify-sklearn \
  --i-classifier ms-classifier.qza \
  --i-reads ms-rep-seqs.qza \
  --o-classification ms-taxonomy.qza

## Next Steps: Taxonomic analysis
