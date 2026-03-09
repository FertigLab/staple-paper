import anndata as ad

# read backed for efficiency
print('reading atlas anndata...')
adata = ad.read_h5ad("data/Human_HMBA_basalganglia_AIT_pre-print.h5ad")
print('done')

# subset to nucleus accumbens
print('subsetting...')
adata_sub = adata[adata.obs["anatomical_region_merged"] == "NAC"].copy()
print('done')

# save
print('saving subset...')
adata_sub.write_h5ad("data/Human_HMBA_basalganglia_AIT_pre-print_NAC.h5ad")
print('all done.')
