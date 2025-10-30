#!/bin/bash
# QIIME2 16S rRNA: running single and paired-end in parrallel
# Dataset 1: barretts_with_oral (paired end)
# Dataset 1: escc_with_control (single end)

set -e

THREADS=64
DADA2_THREADS=32
CLASSIFIER_THREADS=32

# mamba activate qiime2-amplicon-2024.10

mkdir -p qiime2_analysis
cd qiime2_analysis

# Dataset 1
mkdir -p barretts_fastq
for dir in ../barretts_with_oral/SRR*/; do
    sample=$(basename "$dir")
    cp "${dir}/${sample}_1.fastq.gz" "barretts_fastq/${sample}_1.fastq.gz"
    cp "${dir}/${sample}_2.fastq.gz" "barretts_fastq/${sample}_2.fastq.gz"
done

# Dataset 2
mkdir -p escc_fastq
for dir in ../escc_with_control/SRR*/; do
    sample=$(basename "$dir")
    cp "${dir}/${sample}_1.fastq.gz" "escc_fastq/${sample}_1.fastq.gz"
done

echo "Create manifest files..."

# Manifest for Dataset 1
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > barretts_manifest.tsv
for fq in barretts_fastq/*_1.fastq.gz; do
    sample=$(basename "$fq" _1.fastq.gz)
    fwd=$(readlink -f "barretts_fastq/${sample}_1.fastq.gz")
    rev=$(readlink -f "barretts_fastq/${sample}_2.fastq.gz")
    echo -e "${sample}\t${fwd}\t${rev}" >> barretts_manifest.tsv
done

# Manifest for Dataset 2
echo -e "sample-id\tabsolute-filepath" > escc_manifest.tsv
for fq in escc_fastq/*_1.fastq.gz; do
    sample=$(basename "$fq" _1.fastq.gz)
    fwd=$(readlink -f "escc_fastq/${sample}_1.fastq.gz")
    echo -e "${sample}\t${fwd}" >> escc_manifest.tsv
done

echo "Finished reorganising"

###########
# Dataset 1
analyze_barretts() {
    echo "========== Dataset 1: barretts_with_oral (paired-end) =========="

    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path barretts_manifest.tsv \
      --output-path barretts_demux.qza \
      --input-format PairedEndFastqManifestPhred33V2

    qiime demux summarize \
      --i-data barretts_demux.qza \
      --o-visualization barretts_demux.qzv

    qiime dada2 denoise-paired \
      --i-demultiplexed-seqs barretts_demux.qza \
      --p-trim-left-f 0 \
      --p-trim-left-r 0 \
      --p-trunc-len-f 250 \
      --p-trunc-len-r 250 \
      --p-n-threads ${DADA2_THREADS} \
      --o-table barretts_table.qza \
      --o-representative-sequences barretts_rep-seqs.qza \
      --o-denoising-stats barretts_denoising-stats.qza

    qiime metadata tabulate \
      --m-input-file barretts_denoising-stats.qza \
      --o-visualization barretts_denoising-stats.qzv

    qiime feature-table summarize \
      --i-table barretts_table.qza \
      --o-visualization barretts_table.qzv

    qiime feature-table tabulate-seqs \
      --i-data barretts_rep-seqs.qza \
      --o-visualization barretts_rep-seqs.qzv

    qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences barretts_rep-seqs.qza \
      --o-alignment barretts_aligned.qza \
      --o-masked-alignment barretts_masked-aligned.qza \
      --o-tree barretts_unrooted-tree.qza \
      --o-rooted-tree barretts_rooted-tree.qza

    qiime feature-classifier classify-sklearn \
      --i-classifier /data/yzwang/reference/qiime/silva-138-99-nb-classifier.qza \
      --i-reads barretts_rep-seqs.qza \
      --o-classification barretts_taxonomy.qza \
      --p-n-jobs ${CLASSIFIER_THREADS}

    qiime metadata tabulate \
      --m-input-file barretts_taxonomy.qza \
      --o-visualization barretts_taxonomy.qzv

    echo "Finished: Dataset 1"
}

############
# Dataset 2
analyze_escc() {
    echo "========== Dataset 2: escc_with_control (single end) =========="

    qiime tools import \
      --type 'SampleData[SequencesWithQuality]' \
      --input-path escc_manifest.tsv \
      --output-path escc_demux.qza \
      --input-format SingleEndFastqManifestPhred33V2

    qiime demux summarize \
      --i-data escc_demux.qza \
      --o-visualization escc_demux.qzv

    qiime dada2 denoise-single \
      --i-demultiplexed-seqs escc_demux.qza \
      --p-trim-left 0 \
      --p-trunc-len 250 \
      --p-n-threads ${DADA2_THREADS} \
      --o-table escc_table.qza \
      --o-representative-sequences escc_rep-seqs.qza \
      --o-denoising-stats escc_denoising-stats.qza

    qiime metadata tabulate \
      --m-input-file escc_denoising-stats.qza \
      --o-visualization escc_denoising-stats.qzv

    qiime feature-table summarize \
      --i-table escc_table.qza \
      --o-visualization escc_table.qzv

    qiime feature-table tabulate-seqs \
      --i-data escc_rep-seqs.qza \
      --o-visualization escc_rep-seqs.qzv

    qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences escc_rep-seqs.qza \
      --o-alignment escc_aligned.qza \
      --o-masked-alignment escc_masked-aligned.qza \
      --o-tree escc_unrooted-tree.qza \
      --o-rooted-tree escc_rooted-tree.qza

    qiime feature-classifier classify-sklearn \
      --i-classifier /data/yzwang/reference/qiime/silva-138-99-nb-classifier.qza \
      --i-reads escc_rep-seqs.qza \
      --o-classification escc_taxonomy.qza \
      --p-n-jobs ${CLASSIFIER_THREADS}

    qiime metadata tabulate \
      --m-input-file escc_taxonomy.qza \
      --o-visualization escc_taxonomy.qzv

    echo "Finshed: Dataset 2"
}

echo "========== Analysing two datasets in parallel =========="
analyze_barretts &
PID1=$!
analyze_escc &
PID2=$!

wait $PID1
wait $PID2

echo "========== Finshed all datasets =========="

# Merge in different taxonomies
echo "========== Merging taxonomies =========="

mkdir -p merged_results

# Dataset 1
# Family (Level 5)
qiime taxa collapse \
  --i-table barretts_table.qza \
  --i-taxonomy barretts_taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table barretts_table_family.qza

# Genus (Level 6)
qiime taxa collapse \
  --i-table barretts_table.qza \
  --i-taxonomy barretts_taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table barretts_table_genus.qza

# Species (Level 7)
qiime taxa collapse \
  --i-table barretts_table.qza \
  --i-taxonomy barretts_taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table barretts_table_species.qza

# Dataset 2
# Family (Level 5)
qiime taxa collapse \
  --i-table escc_table.qza \
  --i-taxonomy escc_taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table escc_table_family.qza

# Genus (Level 6)
qiime taxa collapse \
  --i-table escc_table.qza \
  --i-taxonomy escc_taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table escc_table_genus.qza

# Species (Level 7)
qiime taxa collapse \
  --i-table escc_table.qza \
  --i-taxonomy escc_taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table escc_table_species.qza

# Combine
echo "Combining..."

# Family (Level 5)
qiime feature-table merge \
  --i-tables barretts_table_family.qza \
  --i-tables escc_table_family.qza \
  --o-merged-table merged_results/merged_table_family.qza

qiime feature-table summarize \
  --i-table merged_results/merged_table_family.qza \
  --o-visualization merged_results/merged_table_family.qzv

# Genus (Level 6)
qiime feature-table merge \
  --i-tables barretts_table_genus.qza \
  --i-tables escc_table_genus.qza \
  --o-merged-table merged_results/merged_table_genus.qza

qiime feature-table summarize \
  --i-table merged_results/merged_table_genus.qza \
  --o-visualization merged_results/merged_table_genus.qzv

# Species (Level 7)
qiime feature-table merge \
  --i-tables barretts_table_species.qza \
  --i-tables escc_table_species.qza \
  --o-merged-table merged_results/merged_table_species.qza

qiime feature-table summarize \
  --i-table merged_results/merged_table_species.qza \
  --o-visualization merged_results/merged_table_species.qzv

#Export merged matrix
echo "========== Exporting merged matrix =========="

# Family (Level 5)
qiime tools export \
  --input-path merged_results/merged_table_family.qza \
  --output-path merged_results/exported_family

biom convert \
  -i merged_results/exported_family/feature-table.biom \
  -o merged_results/feature_table_family.tsv \
  --to-tsv

# Genus (Level 6)
qiime tools export \
  --input-path merged_results/merged_table_genus.qza \
  --output-path merged_results/exported_genus

biom convert \
  -i merged_results/exported_genus/feature-table.biom \
  -o merged_results/feature_table_genus.tsv \
  --to-tsv

# Species (Level 7)
qiime tools export \
  --input-path merged_results/merged_table_species.qza \
  --output-path merged_results/exported_species

biom convert \
  -i merged_results/exported_species/feature-table.biom \
  -o merged_results/feature_table_species.tsv \
  --to-tsv
