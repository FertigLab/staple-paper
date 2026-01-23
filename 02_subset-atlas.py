import anndata as ad

# read backed for efficiency
adata = ad.read_h5ad("data/Human_HMBA_basalganglia_AIT_pre-print.h5ad", backed="r")

# subset to nucleus accumbens
adata = adata[adata.obs["anatomical_region_merged"] == "NAC"]

# save
adata.write_h5ad("data/Human_HMBA_basalganglia_AIT_pre-print_NAC.h5ad")