indir=~/data/metagenomics
outdir=~/metagenomics

mkdir -p $outdir
cd $outdir

conda activate moshpit
qiime info

mosh assembly assemble-megahit --i-reads $indir/reads.qza --p-presets meta-large \
    --p-num-cpu-threads 2 --p-min-contig-len 1000 --o-contigs $outdir/contigs.qza --verbose
    
# N50: A metric that represents the contiguity of a genome assembly. It’s defined as the length of the contig (or scaffold) at which 50% of the entire genome is covered by contigs of that length or longer - the higher this number, the better.
# L50: A metric that represents the number of contigs required to cover 50% of the genome’s total length - the smaller this number, the better.

# assess the quality of contigs assembled in the previous step:
mosh assembly evaluate-quast --i-contigs $outdir/contigs.qza --i-references $indir/reference-genomes.qza \
    --p-threads 4 --p-min-contig 500 --o-reference-genomes $outdir/ref-genomes.qza \
    --o-results-table $outdir/quast-results.qza --o-visualization $outdir/contigs-qc-quast.qzv --verbose
    
#################################################################################

indir=~/data/qiime2
outdir=~/qiime2

mkdir -p $outdir
cd $outdir

conda activate qiime2
qiime metadata tabulate --m-input-file $indir/sample-metadata.tsv --o-visualization sample-metadata.qzv

qiime tools view sample-metadata.qzv
qiime tools extract --input-path sample-metadata.qzv --output-path tmp/

# already imported and de-multiplexed data
qiime demux summarize --i-data $indir/demux.qza --o-visualization demux.qzv
qiime tools view demux.qzv


qiime dada2 denoise-paired --i-demultiplexed-seqs $indir/demux.qza \
  --p-trim-left-f 0 --p-trunc-len-f 250 --p-trim-left-r 0 --p-trunc-len-r 250 \
  --o-representative-sequences asv-seqs.qza --o-table asv-table.qza \
  --o-denoising-stats denoising-stats.qza --o-base-transition-stats base-transition-stats.qza

qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv
qiime tools view denoising-stats.qzv

qiime feature-table summarize \
  --i-table asv-table.qza \
  --m-metadata-file $indir/sample-metadata.tsv \
  --o-summary asv-table.qzv \
  --o-sample-frequencies sample-frequencies.qza \
  --o-feature-frequencies asv-frequencies.qza

qiime feature-table tabulate-seqs \
  --i-data asv-seqs.qza \
  --m-metadata-file asv-frequencies.qza \
  --o-visualization asv-seqs.qzv
qiime tools view asv-seqs.qzv

# filter out any sequences that only appear in one sample (two steps)
# step 1
qiime feature-table filter-features \
  --i-table asv-table.qza \
  --p-min-samples 2 \
  --o-filtered-table asv-table-ms2.qza

# step 2
qiime feature-table filter-seqs \
  --i-data asv-seqs.qza \
  --i-table asv-table-ms2.qza \
  --o-filtered-data asv-seqs-ms2.qza

# train naive Bayues taxonomy classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $indir/reference-sequences.qza \
  --i-reference-taxonomy $indir/reference-taxonomy.qza \
  --o-classifier suboptimal-16S-rRNA-classifier.qza

# apply classifier
qiime feature-classifier classify-sklearn \
  --i-classifier suboptimal-16S-rRNA-classifier.qza \
  --i-reads asv-seqs-ms2.qza \
  --o-classification taxonomy.qza
 
# view classifications
qiime feature-table tabulate-seqs \
  --i-data asv-seqs-ms2.qza \
  --i-taxonomy Greengenes_13_8:taxonomy.qza \
  --m-metadata-file asv-frequencies-ms2.qza \
  --o-visualization asv-seqs-ms2.qzv
qiime tools view asv-seqs-ms2.qzv

# alpha-rarefaction plot]
qiime diversity alpha-rarefaction \
  --i-table asv-table-ms2.qza \
  --p-max-depth 260 \
  --m-metadata-file $indir/sample-metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
qiime tools view alpha-rarefaction.qzv

qiime boots kmer-diversity \
  --i-table asv-table-ms2.qza \
  --i-sequences asv-seqs-ms2.qza \
  --m-metadata-file $indir/sample-metadata.tsv \
  --p-sampling-depth 96 \
  --p-n 10 \
  --p-replacement \
  --p-alpha-average-method median \
  --p-beta-average-method medoid \
  --output-dir boots-kmer-diversity

qiime taxa barplot \
  --i-table asv-table-ms2.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file $indir/sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
qiime tools view  taxa-bar-plots.qzv
 
 
qiime feature-table filter-samples \
  --i-table asv-table-ms2.qza \
  --m-metadata-file $indir/sample-metadata.tsv \
  --p-where '[SampleType] IN ("Human Excrement Compost", "Human Excrement", "Food Compost")' \
  --o-filtered-table asv-table-ms2-dominant-sample-types.qza
  
qiime composition ancombc2 \
  --i-table asv-table-ms2-dominant-sample-types.qza \
  --m-metadata-file $indir/sample-metadata.tsv \
  --p-fixed-effects-formula SampleType \
  --p-reference-levels 'SampleType::Human Excrement Compost' \
  --o-ancombc2-output ancombc2-results.qza
  
qiime composition da-barplot \
  --i-data ancombc2-results.qza \
  --i-taxonomy taxonomy.qza \
  --o-visualization ancombc2-barplot.qzv
qiime tools view  ancombc2-barplot.qzv

  