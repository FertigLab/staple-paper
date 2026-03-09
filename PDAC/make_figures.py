import anndata as ad
import matplotlib
import squidpy as sq
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import json
import pandas as pd
import seaborn as sns
import numpy as np


FIG_WIDTH_MM = 45
FIG_HEIGHT_MM = 45
FIG_WIDTH_IN = FIG_WIDTH_MM / 25.4
FIG_HEIGHT_IN = FIG_HEIGHT_MM / 25.4

plt.rcParams.update({
    'font.size': 5,
    'axes.titlesize': 5,
    'axes.labelsize': 5,
    'xtick.labelsize': 5,
    'ytick.labelsize': 5,
    'legend.fontsize': 5
})

data_dir = 'data/d69909bf-a3c6-48ba-8a5b-2876736de0dc/data/staple'
# Load dataset for sample 4
adata_path = f'{data_dir}/HC04BTC_visiumHD/staple.h5ad'
adata = ad.read_h5ad(adata_path)

#drop NA cells
adata = adata[adata.obs['cell_type']!='NA', :]

# palette in case we need it (uncomment the palette argument in spatial_scatter to use it)
cols = [
    '#f7f7fb',  # ACINAR
    '#f2f0f7',  # B CELLS
    '#eceaf5',  # CYCLING MYELOID
    '#e6e3f2',  # CYCLING TNK
    '#e0ddf0',  # ENDOCRINE
    '#dadaeb',  # ENDOTHELIAL
    '#c7e9c0',  # FIBROBLASTS (light green)
    '#d3d3e8',  # MAST
    '#cdcde6',  # MYELOID
    '#fdae6b',  # PDAC (light orange)
    '#c7c6e3',  # PERICYTES
    '#c1c0e1',  # PLASMA
    '#bbbade',  # TNK
]

palette = mcolors.ListedColormap(cols)

# squidpy spatial scatter
spatial_coords = adata.obsm['spatial']
min_x, min_y = spatial_coords[:, 0].min(), spatial_coords[:, 1].min()
max_x, max_y = spatial_coords[:, 0].max(), spatial_coords[:, 1].max()
sq.pl.spatial_scatter(
    adata,
    color='cell_type',
    size=1.5,
    alpha=0.8,
    #palette=palette,
    palette='tab20c',
    figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
    dpi=300,
    library_id='HC04BTC_visiumHD_hires_image',
    crop_coord=(min_x, min_y, max_x, max_y),
    legend_loc='right margin'
)

spatial_fig = plt.gcf()
spatial_ax = plt.gca()
legend = spatial_ax.get_legend()

if legend is not None:
    handles, labels = spatial_ax.get_legend_handles_labels()
    legend.remove()
    spatial_fig.savefig('HC04BTC_spatial_scatter.png', dpi=300, bbox_inches='tight')

    legend_fig = plt.figure(figsize=(FIG_WIDTH_IN * 0.5, FIG_HEIGHT_IN), dpi=300)
    legend_ax = legend_fig.add_subplot(111)
    legend_ax.axis('off')
    legend_ax.legend(handles, labels, loc='center left', frameon=False, fontsize=5)
    legend_fig.savefig('HC04BTC_spatial_scatter_legend.svg', bbox_inches='tight')
    plt.close(legend_fig)
else:
    spatial_fig.savefig('HC04BTC_spatial_scatter.png', dpi=300, bbox_inches='tight')

plt.close(spatial_fig)


# sc.gr interactions
sq.gr.interaction_matrix(adata, cluster_key="cell_type", normalized=True)
sq.pl.interaction_matrix(adata, cluster_key="cell_type", 
                         map='viridis', 
                         figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
                         dpi=300,
                         save='HC04BTC_interaction_matrix.png'
                         )


# spacemarkers LR scores
df = adata.uns['LRscores'].loc[:,['FIBROBLASTS_to_PDAC','FIBROBLASTS_to_ACINAR']].sort_values(by='FIBROBLASTS_to_PDAC', ascending=False).head(5)
df.insert(1, '...', '...')
df.to_csv('LR_scores.csv')


# xsample spatial neighborhoods
nbr = json.load(open(f"{data_dir}/reports/mqc/neighbors_mqc.json"))
nbr = pd.DataFrame(nbr['data'][1])
nbr = nbr[nbr.index!='NA']

# horizontal stacked barplot of nbr using the same cell_type colors as spatial_scatter
nbr_plot = nbr.sort_index().transpose().sort_index(ascending=False)

# Ensure colors are matched by cell type name, not by positional order.
cell_type_order = list(adata.obs['cell_type'].astype('category').cat.categories)
cell_type_color_map = dict(zip(cell_type_order, adata.uns['cell_type_colors']))
nbr_colors = [cell_type_color_map.get(col, '#bdbdbd') for col in nbr_plot.columns]
nbr_plot.index = nbr_plot.index.str.replace('_visiumHD', '')  # Clean up index labels

ax = nbr_plot.plot(
    kind='barh',
    stacked=True,
    color=nbr_colors,
    figsize=(FIG_WIDTH_IN*1.65, FIG_HEIGHT_IN),
    legend=False,
    xlabel='Number of Neighbors',
    width=0.95
)
ax.set_xticks(np.arange(0, 500000, 100000))
ax.set_title('FIBROBLASTS Neighbors by Cell Type', fontweight='bold')
plt.tight_layout()
plt.savefig('HC04BTC_spatial_neighbors.svg', dpi=300, bbox_inches='tight')
plt.close()

# xsample cell spacemarkers ligrec scores
ligrec = json.load(open(f"{data_dir}/reports/mqc/lrscores_diff_response_results_mqc.json"))
ligrec = pd.DataFrame(ligrec['data'])
ligrec.index = ligrec.index.str.replace('-FIBROBLASTS_to_PDAC', '')  # Clean up index labels
ligrec = ligrec.transpose().sort_index()
ligrec.index = ligrec.index.str.replace('_visiumHD', '')  # Clean up index labels

#heatmap
plt.figure(figsize=(FIG_WIDTH_IN*2, FIG_HEIGHT_IN*1.25), dpi=300)
cmap = nbr_colors[4:8]  # Use a subset of the same colors for the heatmap
cmap.reverse()  # Reverse the color map to have higher scores in darker colors
sns.heatmap(ligrec, cmap=cmap, cbar_kws={'label': 'Interaction Score'})
plt.title('FIBROBLASTS to PDAC Interactions', fontweight='bold')
plt.xlabel('Ligand-Receptor Pairs')
plt.tight_layout()
#make colors pastel
plt.savefig('HC04BTC_ligrec_heatmap.svg', dpi=300, bbox_inches='tight')
plt.close() 