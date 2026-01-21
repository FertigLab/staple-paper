
# download the raw archive from GSE307586
import os
import urllib.request
import tarfile
import pandas as pd
import numpy as np
import gzip
import shutil
import h5py
import scanpy as sc



# Define the GSE accession number and construct the download URL
gse = "GSE307586"
url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse}&format=file"

# Create the output directory if it doesn't exist
output_dir = f"data/{gse}"
os.makedirs(output_dir, exist_ok=True)

# Download the file if it doesn't already exist
response = urllib.request.urlopen(url)
response_url = response.geturl()
output_file = os.path.join(output_dir, f"{gse}.tar")
urllib.request.urlretrieve(response_url, output_file) if not os.path.exists(output_file) else None

# Extract the downloaded tar file
with tarfile.open(output_file, "r") as tar:
    tar.extractall(path=f"{output_dir}/extracted")

# Parse individual filenames and put them into sample folders unzipped
all_files = os.listdir(f"{output_dir}/extracted")
samples = np.unique([x.split("_")[0] for x in all_files])

# track theses names and place them in the sample folders
raw_names = ["tissue_hires_image.png",
               "tissue_lowres_image.png",
               "aligned_fiducials.jpg",
               "detected_tissue_image.jpg",
               "scalefactors_json.json",
               "tissue_positions_list.csv",
               "barcodes.tsv.gz",
               "features.tsv.gz",
               "matrix.mtx.gz"]


for s in samples:
    os.makedirs(f"{output_dir}/samples/{s}/filtered_feature_bc_matrix", exist_ok=True)
    os.makedirs(f"{output_dir}/samples/{s}/spatial", exist_ok=True)
    for f in all_files:
        if f.startswith(s):
            #find proper name
            f2 = raw_names[np.where([x in f for x in raw_names])[0][0]]
            #unzip and move to appropriate folders
            if f2 in ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]:
                shutil.copyfile(f"{output_dir}/extracted/{f}", f"{output_dir}/samples/{s}/filtered_feature_bc_matrix/{f2}")
            else:
                with gzip.open(f"{output_dir}/extracted/{f}", 'rb') as f_in:
                    with open(f"{output_dir}/samples/{s}/spatial/{f2}", 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)


# use scanpy to create the 10X filtered_feature_matrix.h5 format 
# files form the matrix.mtx, barcodes.tsv, features.tsv
# gemini provided a function to write 10x h5 files
def write_10x_h5(adata, filename):
    # 1. 10x format requires CSC (Compressed Sparse Column)
    # and the matrix should be (features x barcodes)
    # We transpose the AnnData (cells x features) to get (features x cells)
    X_csc = adata.X.T.tocsc()

    with h5py.File(filename, 'w') as f:
        group = f.create_group("matrix")
        
        # Datasets expected by 10x readers
        group.create_dataset("barcodes", data=adata.obs_names.values.astype('S'))
        group.create_dataset("data", data=X_csc.data)
        group.create_dataset("indices", data=X_csc.indices)
        group.create_dataset("indptr", data=X_csc.indptr)
        
        # 10x format stores shape as [features, barcodes]
        group.create_dataset("shape", data=np.array(X_csc.shape).astype(np.int32))
        
        # Feature metadata group
        features = group.create_group("features")
        features.create_dataset("id", data=adata.var_names.values.astype('S'))
        features.create_dataset("name", data=adata.var_names.values.astype('S'))
        features.create_dataset("feature_type", data=np.array(["Gene Expression"] * adata.n_vars).astype('S'))
        features.create_dataset("genome", data=np.array(["GRCh38"] * adata.n_vars).astype('S'))


for s in samples:
    adata = sc.read_10x_mtx(f"{output_dir}/samples/{s}/filtered_feature_bc_matrix", var_names='gene_symbols', cache=False)
    write_10x_h5(adata, f"{output_dir}/samples/{s}/filtered_feature_bc_matrix.h5")

#test
from spatialdata_io import visium
from spatialdata_io.experimental import to_legacy_anndata
for s in samples:
    try:
        sdata = visium(f"{output_dir}/samples/{s}", dataset_id=s)
        adata = to_legacy_anndata(sdata, coordinate_system=s)
        print(adata)
    except Exception as e:
        print(f"Error loading sample {s}: {e}")
print("All done.")
