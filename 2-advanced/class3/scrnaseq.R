# title: "Seurat seurat tutorial"
# output: html_document
# date: "2026-10-19"
# URL: "https://satijalab.org/seurat/articles/pbmc3k_tutorial"
# Essential commands: "https://satijalab.org/seurat/articles/essential_commands"

outdir <- '~/scrnaseq'
if (!dir.exists(outdir))
  dir.create(outdir, recursive = TRUE)
knitr::opts_knit$set(root.dir = outdir)
knitr::opts_chunk$set(echo = TRUE)

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Load the seurat dataset
# Seuratのデータセットを読み込みます。
pbmc_data <- Read10X(data.dir = "~/data/scrnaseq/filtered_gene_bc_matrices/hg19/")

# Lets examine a few genes (rows) and cells (columns)
# いくつかの遺伝子（行）と細胞（列）を確認します。
pbmc_data[c("CD3D", "TCL1A", "MS4A1"), c("AAACATACAACCAC-1", "AAACATTGAGCTAC-1", "AAACATTGATCAGC-1")]

# Initialize the Seurat object with the raw (non-normalized) count data.
# 生の（正規化されていない）カウントデータを使ってSeuratオブジェクトを初期化します。

# Note that we can set some baseline filters to exclude genes expressed in less than 3 cells and cells with less than 200 genes
# 基本的なフィルタを設定して、3個未満の細胞でしか発現していない遺伝子と、検出遺伝子数が200未満の細胞を除外できます。

seurat <- CreateSeuratObject(counts = pbmc_data, project = "pbmc3k", min.cells = 3, min.features = 200)
seurat

str(seurat)

# This type of R object stores information in slots using the @ symbol, such as seurat@assays and seurat@metadata.
# この種類のRオブジェクトでは、seurat@assaysやseurat@metadataのように、@記号を使ったスロットに情報が保存されます。

# However, most data is readily accessible using shortcuts using the [["name"]] notation
# ただし、多くのデータには[["name"]]という記法を使って簡単にアクセスできます。
head(seurat[["RNA"]]$counts)

# Low-quality / dying cells often exhibit extensive mitochondrial contamination
# 品質の低い細胞や死にかけた細胞では、ミトコンドリア由来のシグナルが多く見られることがあります。

# We can add our own columns, such as percent.mt, as metadata using the [["percent.mt"]] notation
# [["percent.mt"]]という記法を使い、percent.mtなどの独自の列をメタデータとして追加できます。
seurat[["percent_MT"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# Show QC metrics for the first 5 cells
# 最初の5細胞の品質管理（QC）指標を表示します。
head(seurat@meta.data, 5)

# Visualize QC metrics as a violin plot
# QC指標をバイオリンプロットで可視化します。
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), layer='counts', ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# FeatureScatterは通常、特徴量同士の関係を可視化するために使いますが、次のような用途にも使えます。
# オブジェクトで計算された値、たとえばメタデータの列や主成分（PC）スコアなどです。
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent_MT")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# We filter cells that have unique feature counts over 2,500 or less than 200 and >5% mitochondrial counts
# 検出された特徴量数が2,500を超える細胞、200未満の細胞、またはミトコンドリア比率が5%を超える細胞を除外します。
print(paste0('Before filtering: ', ncol(seurat[["RNA"]]$counts), ' features'))

# In addition to built-in fields such as nFeature_RNA, we can also use our custom metadata columns in the subset() function
# nFeature_RNAなどの組み込みフィールドに加え、subset()関数では独自に追加したメタデータ列も利用できます。
seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_MT < 5)
print(paste0('After filtering: ', ncol(seurat[["RNA"]]$counts), ' features'))

# We normalize the feature expression measurements for each cell by the total expression,
# multiply by a scale factor (10,000 by default), and log-transforms the result.
# 各細胞の特徴量発現量を、その細胞の総発現量で正規化し、スケーリング係数（既定値は10,000）を掛けてから、対数変換します。
# See SCTransform() for a more sophisticated normalization method
# より高度な正規化方法については、SCTransform()を参照してください。
#seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- NormalizeData(seurat)

# Our analysis should mainly focus on genes that vary in expression levels among cells.
# 解析では、細胞間で発現量が変動する遺伝子を主に対象とします。
# Calculate a subset of genes that exhibit high cell-to-cell variation in the dataset.
# データセット内で細胞間変動が大きい遺伝子の一部を選びます。
# These genes are highly expressed in some cells but are have low expression in others.
# これらの遺伝子は、ある細胞では高く発現し、別の細胞では低く発現します。
# The top 2000 most variable genes will be selected and used for downstream analysis
# 変動性が最も高い上位2,000遺伝子を選び、以降の解析に使用します。
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000, verbose=FALSE)

# List the 10 most highly variable genes
# 変動性が最も高い上位10遺伝子を表示します。
top10 <- head(VariableFeatures(seurat), 10)
top10

# Plot variable genes along with their labels (dispersion vs mean plot)
# 変動遺伝子をラベル付きでプロットします（分散と平均のプロット）。
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Scale the expression of each gene, so that the mean expression across cells is 0 and the variance is 1,
# so that highly-expressed genes do not dominate. The results of this are stored in seurat[["RNA"]]$scale.data.
# 各遺伝子の発現量をスケーリングし、細胞間での平均を0、分散を1にします。
# これにより発現量の高い遺伝子だけが解析を支配するのを防ぎます。結果はseurat[["RNA"]]$scale.dataに保存されます。

# By default, only variable features are scaled, but you can specify the features argument to scale additional features.
# 既定では変動特徴量のみがスケーリングされますが、features引数で対象を追加できます。
seurat <- ScaleData(seurat, features = rownames(seurat))

# The RNA assay now contains three layers, count, data, and scale.data
# RNAアッセイには、counts、data、scale.dataの3つのレイヤーが含まれています。
print('Raw counts')
seurat[["RNA"]]$counts[5:10,]

print('Normalized counts')
seurat[["RNA"]]$data[5:10,]

print('Scaled normalized counts')
seurat[["RNA"]]$scale.data[5:10, 1:10]

# Perform linear dimensional reduction using PCA. Note that only the highly variable genes we identified above are used for PCA
# PCAを使って線形次元削減を行います。PCAには、先ほど特定した高変動遺伝子のみを使用します。
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))

# PCA data is stored in the seurat@reductions$pca slot but can be accessed easily through seurat[["pca"]]
# PCAのデータはseurat@reductions$pcaスロットに保存されますが、seurat[["pca"]]から簡単にアクセスできます。
#print(seurat[["pca"]], dims = 1:5, nfeatures = 5)

# Let's examine and visualize PCA results a few different ways.
# PCAの結果をさまざまな方法で確認・可視化します。

# Plot the first two principal components (dims 1 and 2). Which genes contribute most?
# 最初の2つの主成分（dims 1と2）をプロットします。どの遺伝子の寄与が大きいでしょうか？
VizDimLoadings(seurat, dims = 1:2, reduction = "pca")

# A limitation of PCA is that we can only plot two components at a time. This plot displays components 1 vs 2.
# PCAでは一度に2つの成分しかプロットできません。この図では第1主成分と第2主成分を表示します。
DimPlot(seurat, reduction = "pca") + NoLegend()

# A heatmap clearly shows the expression patterns of the most influential genes within each listed principal component (1, 2, and 3)
# ヒートマップでは、指定した各主成分（1、2、3）に強く寄与する遺伝子の発現パターンを確認できます。
DimHeatmap(seurat, dims = c(1,2,3), cells = 500, balanced = TRUE)

# PCA reduces the number of dimensions, but we still need to decide how many components to include.
# PCAで次元数を減らした後、いくつの主成分を解析に含めるか決める必要があります。
# The "elbow" in the plot below indicates a point of diminishing returns at around 10, which is a reasonable number.
# 下のプロットの「肘（elbow）」は、効果の増加が小さくなり始める点を示します。ここでは約10個が妥当な目安です。
# However, you may want to re-run the analysis using a different number of components
# ただし、主成分数を変えて解析を再実行してもよいでしょう。
ElbowPlot(seurat)

# We first construct a KNN graph based on the Euclidean distance in PCA space using the number of components determined above (10).
# まず、上で決めた10個の主成分を使い、PCA空間でのユークリッド距離に基づくK近傍（KNN）グラフを作成します。
# The edge weights between any two cells are refined based on the shared overlap in their local neighborhoods.
# 細胞間のエッジの重みは、それぞれの局所近傍がどの程度重複するかに基づいて調整されます。
num_principal_components <- 10
seurat <- FindNeighbors(seurat, dims = 1:num_principal_components)

# Next the cells are iteratively clustered together using a optimization techniques such as the Louvain or Leiden algorithms.
# 次に、Louvain法やLeiden法などの最適化アルゴリズムを使って、細胞を反復的にクラスタリングします。
# The value of the resolution parameter affects the number of clusters and should be adjusted if necessary (0.4-1.2)
# resolutionパラメータはクラスタ数に影響するため、必要に応じて調整します（目安：0.4～1.2）。
resolution <- 0.5
seurat <- FindClusters(seurat, resolution = resolution)

# Look at cluster IDs of the first 5 cells
# 最初の5細胞のクラスタIDを確認します。
head(Idents(seurat), 5)

# While PCA reduced the number of dimensions, it was hard to understand tha structure since we could only plot two components at a time.
# PCAで次元数を減らしても、一度に2成分しか表示できないため、データ構造を理解するのは容易ではありません。

# UMAP and t-SNE are 2D visualization techniques that help ensure that similar cells are displayed together.
# UMAPとt-SNEは、似た細胞が近くに表示されるようにする2次元可視化手法です。

# This function calculates UMAP coordinates but does not plot them.
# この関数はUMAP座標を計算しますが、プロット自体は作成しません。
seurat <- RunUMAP(seurat, dims = 1:num_principal_components)

# Plot the calculated UMAP coordinates (stored in seurat@reductions[["umap"]])
# 計算したUMAP座標（seurat@reductions[["umap"]]に保存）をプロットします。
DimPlot(seurat, reduction = "umap", label=TRUE)

# At this point, we can save the workspace session as an rds file and reload it later or share it with a collaborator
# ここで作業環境をRDSファイルとして保存できます。後で読み込んだり、共同研究者と共有したりできます。
rds_file <- paste0(outdir, '/pbmc_tutorial.rds')
saveRDS(seurat, file = rds_file)

# Later, you can reload the the workspace session using the readRDS function
# 後でreadRDS関数を使って作業環境を再読み込みできます。
# seurat <- readRDS(rds_file)

###############################################################
# Find cluster-specific markers
# クラスター特異的なマーカーを特定する
###############################################################

# Each cluster is generally assumed to represent a different cell type or sub-type with distinct marker genes.
# 通常、各クラスターは、特徴的なマーカー遺伝子を持つ異なる細胞型またはサブタイプに対応すると考えます。
# FindMarkers tries to find unique markers that separate the selected cluster from all others.
# FindMarkersは、選択したクラスターを他のすべてのクラスターから区別するマーカーを探索します。
# Note: This takes a few minutes to run!
# 注意：実行には数分かかります。
cluster3.markers <- FindMarkers(seurat, ident.1 = 3)

# Display the list of marker genes
# マーカー遺伝子の一覧を表示します。
head(cluster3.markers, n = 5)

# However, it is often useful to find markers that distinguish specific clusters from each other (cluster 5 vs clusters 0 and 3)
# 特定のクラスター同士を区別するマーカーを調べることも有用です（例：クラスター5とクラスター0・3の比較）。
FindMarkers(seurat, ident.1 = 5, ident.2 = c(0, 3))

# You can plot expression levels of specific genes across all clusters
# 特定の遺伝子の発現量を、すべてのクラスターについてプロットできます。
VlnPlot(seurat, features = c("PPBP", "FCER1A", "CD79A"))

# FeaturePlot plots the candidate markers on the UMAP plot instead
# FeaturePlotを使うと、候補マーカーの発現をUMAP上に表示できます。
FeaturePlot(seurat, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# Find markers for every cluster compared to all remaining cells, report only the positive ones
# 各クラスターを残りすべての細胞と比較してマーカーを探索し、正のマーカーのみを報告します。
# Takes a long time!
# 実行に時間がかかります！
pbmc.markers <- FindAllMarkers(seurat, only.pos = TRUE)

# Group the markers by cluster and keep only markers with log2 fold change greater than 1
# マーカーをクラスターごとにまとめ、log2 fold changeが1を超えるものだけを残します。
pbmc.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)

# However, a low p value does not necessarily make it a good marker
# ただし、p値が小さいからといって、必ずしも良いマーカーとは限りません。
VlnPlot(seurat, features = c("LDHB"))

# However, an expression heatmap of the top 10 markers for each cluster shows how each cluster has a distinct epxression pattern
# 各クラスターの上位10個のマーカー遺伝子の発現ヒートマップを見ると、クラスターごとに異なる発現パターンが分かります。
pbmc.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup() -> top10
DoHeatmap(seurat, features = top10$gene) + NoLegend()

###############################################################
# Automatic cell type annotation with SingleR
# SingleRによる細胞型の自動アノテーション
###############################################################

library(SingleR)
library(celldex)
library(SingleCellExperiment)

# Make sure the clusters are stored explicitly (Idents may get overwritten later)
# クラスター情報を明示的に保存しておきます（後でIdentsが上書きされる可能性があります）。
seurat$cluster <- as.character(seurat$seurat_clusters)

# Load reference datasets from celldex ---------------------------------
# celldexから参照データセットを読み込みます。
# - HumanPrimaryCellAtlasData: broad cell types (e.g. "T_cells", "Monocyte")
# - HumanPrimaryCellAtlasData：大まかな細胞型（例："T_cells"、"Monocyte"）
# - MonacoImmuneData: fine-grained immune types (e.g. "Naive CD4 T cells")
# - MonacoImmuneData：より詳細な免疫細胞型（例："Naive CD4 T cells"）
ref.hpca   <- celldex::HumanPrimaryCellAtlasData()
ref.monaco <- celldex::MonacoImmuneData()

# Inspect what labels are available
# 利用可能なラベルを確認します。
print('Human Primary Cell Atlas')
table(ref.hpca$label.main)

print('Monaco Immune Data coarse cell types')
table(ref.monaco$label.main)

print('Monaco Immune Data fine cell types')
table(ref.monaco$label.fine)

# Convert Seurat object -> SingleCellExperiment -----------------------
# SeuratオブジェクトをSingleCellExperimentに変換します。
# SingleR expects a SingleCellExperiment object with log-normalized expression data ('logcounts')
# SingleRは、対数正規化された発現データ（'logcounts'）を含むSingleCellExperimentオブジェクトを必要とします。
DefaultAssay(seurat) <- "RNA"
pbmc.sce <- as.SingleCellExperiment(seurat)
pbmc.sce  # check that "logcounts" is present in assays()

# Single-cell level annotation: Each cell is scored against the reference independently.
# 細胞単位のアノテーション：各細胞を参照データと個別に照合してスコアを計算します。
pred.hpca.cell <- SingleR(test = pbmc.sce, ref = ref.hpca, labels = ref.hpca$label.main)
table(pred.hpca.cell$labels)

# Now annotate each cell using Monaco Fine reference data
# 次に、Monacoの詳細な参照データを使って各細胞をアノテーションします。
pred.monaco.cell <- SingleR(test = pbmc.sce, ref = ref.monaco,  labels = ref.monaco$label.fine)
table(pred.monaco.cell$labels)

# Apply a label to each cluster instead of cell using HPCA
# HPCAを使い、細胞ごとではなくクラスターごとにラベルを付けます。
pred.hpca.clust <- SingleR(test = pbmc.sce, ref = ref.hpca, labels = ref.hpca$label.main, clusters = seurat$cluster)
pred.hpca.clust$labels

# Apply a label to each cluster instead of cell using Monaco fine
# Monacoの詳細ラベルを使い、細胞ごとではなくクラスターごとにラベルを付けます。
pred.monaco.clust <- SingleR(test = pbmc.sce, ref = ref.monaco, labels = ref.monaco$label.fine, clusters = seurat$cluster)
pred.monaco.clust$labels

# Add results to the Seurat object
# 結果をSeuratオブジェクトに追加します。
seurat$SingleR_HPCA_cell   <- pred.hpca.cell$labels
seurat$SingleR_HPCA_pruned <- pred.hpca.cell$pruned.labels
seurat$SingleR_Monaco_cell <- pred.monaco.cell$labels

# Per-cluster labels: map each cluster ID to its label
# クラスター単位のラベル：各クラスターIDを対応するラベルに変換します。
hpca.clust.map <- setNames(pred.hpca.clust$labels, rownames(pred.hpca.clust))
monaco.clust.map <- setNames(pred.monaco.clust$labels, rownames(pred.monaco.clust))

print('HCPA cluster names')
print(hpca.clust.map)

print('Monaco cluster names')
print(monaco.clust.map)

seurat$SingleR_HPCA_cluster   <- unname(hpca.clust.map[seurat$cluster])
seurat$SingleR_Monaco_cluster <- unname(monaco.clust.map[seurat$cluster])

# Visualize HPCA cluster names on the UMAP plot
# UMAPプロット上にHPCAのクラスター名を表示します。
p.hpca.cell <- DimPlot(seurat, reduction = "umap", group.by = "SingleR_HPCA_cell", label = TRUE, repel = TRUE) +
  ggtitle("SingleR: HPCA (per cell)") + NoLegend()

p.hpca.clust <- DimPlot(seurat, reduction = "umap", group.by = "SingleR_HPCA_cluster", label = TRUE, repel = TRUE) +
  ggtitle("SingleR: HPCA (per cluster)") + NoLegend()

p.hpca.cell + p.hpca.clust

# Visualize Monaco on the UMAP plot
# UMAPプロット上にMonacoのラベルを表示します。
p.monaco.cell <- DimPlot(seurat, reduction = "umap", group.by = "SingleR_Monaco_cell", label = TRUE, repel = TRUE) +
  ggtitle("SingleR: Monaco (per cell)") + NoLegend()

p.monaco.clust <- DimPlot(seurat, reduction = "umap", group.by = "SingleR_Monaco_cluster", label = TRUE, repel = TRUE) +
  ggtitle("SingleR: Monaco (per cluster)") + NoLegend()

p.monaco.cell + p.monaco.clust
