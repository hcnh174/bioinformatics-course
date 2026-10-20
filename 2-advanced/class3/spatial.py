# 1. Import libraries and set aliases
# ライブラリをインポートし、エイリアスを設定する
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
import anndata as ad
import warnings
from pathlib import Path

# Create the output directory within your home folder 
# ホームフォルダー内に、出力ディレクトリーを作成してください。
outdir = Path.home() / "spatial"
outdir.mkdir(parents=True, exist_ok=True)
os.chdir(outdir)
print(os.getcwd())

# Set default scanpy settings
# Scanpyのデフォルト設定を行う
sc.settings.verbosity = 2
sc.settings.figdir = "figures"
sc.set_figure_params(dpi=100, facecolor="white", figsize=(5, 5))

# 2. Load the data (takes some time to download)
# 2. データを読み込む（ダウンロードに時間がかかる場合があります）
warnings.filterwarnings('ignore')
adata = sq.datasets.visium(sample_id="V1_Breast_Cancer_Block_A_Section_1")
adata.var_names_make_unique()
print(adata)  # AnnData: spots x genes; coordinates in adata.obsm["spatial"]

# Examine the tissue image with spots overlaid but without expression data
# 発現データを含めず、スポットを重ね合わせた組織画像を確認する
sc.pl.spatial(adata, img_key="hires", alpha=0.5, title="Visium spots on H&E")

# 3. Quality control
# 3. 品質管理
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# QC distributions
# QC分布
fig, axs = plt.subplots(1, 3, figsize=(14, 3.5))
sns_kw = dict(bins=50, color="steelblue")
axs[0].hist(adata.obs["total_counts"], **sns_kw); axs[0].set_title("Total counts per spot")
axs[1].hist(adata.obs["n_genes_by_counts"], **sns_kw); axs[1].set_title("Genes per spot")
axs[2].hist(adata.obs["pct_counts_mt"], **sns_kw); axs[2].set_title("% mitochondrial")
plt.tight_layout(); plt.show()

# Spatial QC. Low-count areas may represent real biology (necrosis, fat, acellular stroma) rather than technical error.
# 空間的品質管理（Spatial QC）。カウント数が低い領域は、技術的なエラーではなく、実際の生物学的状態（壊死、脂肪、細胞成分を欠く間質など）を反映している可能性がある。
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts", "pct_counts_mt"])

# 4. Filtering using permissive thresholds
# 4. 許容的な閾値を用いたフィルタリング
print("Spots before:", adata.n_obs)
sc.pp.filter_cells(adata, min_counts=500)
adata = adata[adata.obs["pct_counts_mt"] < 25].copy()
sc.pp.filter_genes(adata, min_cells=10)
print("Spots after:", adata.n_obs, "| genes after:", adata.n_vars)

# 5. Normalize and log-transform counts. Identify highly variable genes
# 5. カウントの正規化および対数変換を行う。変動の大きい遺伝子を特定する。
adata.layers["counts"] = adata.X.copy()          # keep raw counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# 6. Dimensionality reduction using PCA. Identify nearest neighbors and iteratively create clusters
# 6. PCAを用いた次元削減。近傍点を特定し、反復的にクラスタを生成する
sc.pp.pca(adata, n_comps=30)  # uses the highly variable genes by default
sc.pp.neighbors(adata, n_pcs=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.6, key_added="clusters", flavor="igraph", n_iterations=2, directed=False)

# Compare the clusters in UMAP space vs. tissue space
# UMAP空間と組織空間におけるクラスターを比較する
sc.pl.umap(adata, color="clusters", title="Clusters (expression space)")
sc.pl.spatial(adata, img_key="hires", color="clusters", title="Clusters (tissue space)")

# 7. Identify marker genes that distinguish each cluster from all other clusters
# 7. 各クラスターを他のすべてのクラスターから区別するマーカー遺伝子を特定する
sc.tl.rank_genes_groups(adata, groupby="clusters", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata, n_genes=4, standard_scale="var")

# Find and filter markers based on adjusted p-values and log fold change
# 調整済みp値および対数倍率変化（log fold change）に基づいてマーカーを検索・抽出する
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers = markers[(markers["pvals_adj"] < 0.01) & (markers["logfoldchanges"] > 1)]
markers.to_csv("cluster_markers.csv", index=False)
markers.groupby("group").head(3)

# 8. Annotate spatial regions using known marker signatures associated with breast cancer
# 8. 乳がんに伴う既知のマーカーシグネチャーを用いて、空間領域に注釈を付与する
signatures = {
    "Epithelial_tumor": ["EPCAM", "KRT8", "KRT18", "KRT19", "ERBB2"],
    "Proliferation":    ["MKI67", "TOP2A", "PCNA", "MCM2"],
    "Stroma_CAF":       ["COL1A1", "COL1A2", "DCN", "LUM", "FAP"],
    "T_cell":           ["CD3D", "CD3E", "CD2", "IL7R"],
    "B_cell_plasma":    ["MS4A1", "CD79A", "IGKC", "JCHAIN"],
    "Myeloid":          ["CD68", "CD14", "LYZ", "CD163"],
    "Hypoxia":          ["VEGFA", "CA9", "SLC2A1", "LDHA", "BNIP3"],
}

# Score genes with respect to these signatures
# これらのシグネチャーに基づいて遺伝子をスコアリングする
for name, genes in signatures.items():
    present = [g for g in genes if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=present, score_name=name)

# Plot using a distinct color for each signature
# 各シグネチャーに異なる色を使用してプロットする
sc.pl.spatial(adata, img_key="hires", color=list(signatures.keys()), ncols=4, cmap="magma", wspace=0.1)

# Label clusters using the mean signature score
# 平均シグネチャースコアを用いてクラスターにラベルを付ける
score_df = adata.obs.groupby("clusters")[list(signatures.keys())].mean()
fig, ax = plt.subplots(figsize=(7, 4))
im = ax.imshow(score_df.values, cmap="viridis", aspect="auto")
ax.set_xticks(range(score_df.shape[1])); ax.set_xticklabels(score_df.columns, rotation=45, ha="right")
ax.set_yticks(range(score_df.shape[0])); ax.set_yticklabels(score_df.index)
ax.set_ylabel("Cluster"); plt.colorbar(im, label="Mean score"); plt.tight_layout(); plt.show()

# Plot with respect to the first two clusters (0 and 1)
# 最初の2つのクラスター（0および1）に関するプロット
group_a, group_b = "0", "1"
sc.tl.rank_genes_groups(adata, groupby="clusters", groups=[group_a], reference=group_b, method="wilcoxon", key_added="de_a_vs_b")
de = sc.get.rank_genes_groups_df(adata, group=group_a, key="de_a_vs_b")
print(de.head(15))
sc.pl.rank_genes_groups(adata, key="de_a_vs_b", n_genes=15, sharey=False)

# Simple volcano plot comparing clusters 0 and 1
# クラスター0と1を比較するシンプルなボルケーノプロット
de["neglog10_padj"] = -np.log10(de["pvals_adj"].clip(lower=1e-300))
plt.figure(figsize=(5, 4))
plt.scatter(de["logfoldchanges"], de["neglog10_padj"], s=4, c="grey")
sig = de[(de["pvals_adj"] < 0.01) & (de["logfoldchanges"].abs() > 1)]
plt.scatter(sig["logfoldchanges"], sig["neglog10_padj"], s=6, c="crimson")
plt.xlabel(f"log2FC (cluster {group_a} vs {group_b})"); plt.ylabel("-log10 adj. p")
plt.tight_layout(); plt.show()

# 10. Spatially aware analysis with Squidpy
# 10. Squidpyを用いた空間情報を考慮した解析
# 10a. Build the spatial neighbor graph (Visium = hexagonal grid, 6 neighbors)
# 10a. 空間近傍グラフの構築（Visiumの場合：六角形グリッド、近傍数6）
sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)

# 10b. Neighborhood enrichment: which clusters are adjacent more often than expected by chance? (e.g., tumor/immune borders)
# 10b. 近傍の濃縮（Neighborhood enrichment）：偶然期待される頻度よりも高い頻度で隣接しているクラスターはどれか？（例：腫瘍と免疫細胞の境界領域など）
sq.gr.nhood_enrichment(adata, cluster_key="clusters")
sq.pl.nhood_enrichment(adata, cluster_key="clusters", method="average", figsize=(5, 5))

# 10c. Identify spatially variable genes using Moran's I (test for spatial autocorrelation)
# 10c. Moran's I（空間自己相関の検定）を用いて、空間的に変動する遺伝子を特定する
adata_hvg = adata[:, adata.var["highly_variable"]].copy()   # speed: HVGs only
sq.gr.spatial_autocorr(adata_hvg, mode="moran", n_perms=100, n_jobs=1, seed=0)
svg = adata_hvg.uns["moranI"].sort_values("I", ascending=False)
print(svg.head(15))
svg.to_csv("spatially_variable_genes.csv")

# Examine the top 6 identified spatially variable genes
# 特定された空間変動遺伝子のうち、上位6つを検討する
top_svg = svg.head(6).index.tolist()
sc.pl.spatial(adata, img_key="hires", color=top_svg, ncols=3, cmap="viridis")

# 11. Save the adata object in h5ad format
# 11. adataオブジェクトをh5ad形式で保存する
adata.write_h5ad("breast_cancer_visium_processed.h5ad")
print("Done. Outputs: processed .h5ad, cluster_markers.csv, spatially_variable_genes.csv")

# 12. Later, the h5ad file can be reloaded
# 12. 後から h5ad ファイルを再読み込みできます
#adata = ad.read_h5ad("breast_cancer_visium_processed.h5ad")

