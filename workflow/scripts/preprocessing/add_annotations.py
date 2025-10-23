# %% 
import pandas as pd
import scanpy as sc
import os
import platform

# %% Parse snakemake
if 'snakemake' in locals():
    h5ad_input = snakemake.input['h5ad_input']
    annotation_dict = snakemake.input['annotation_dict']
    h5ad_output = snakemake.output['h5ad_output']
else:
    h5ad_input = 'results/preprocessing/Julio_OCEAN_Nereid.h5ad'
    annotation_dict = 'data/OCEAN_LR.csv'
    h5ad_output = 'results/preprocessing/Julio_OCEAN_Nereid.h5ad'
    if platform.system() == 'Linux':
        h5ad_input = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', h5ad_input)
        h5ad_output = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', h5ad_output)
        annotation_dict = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', annotation_dict)

    elif platform.system() == 'Darwin':
        h5ad_input = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', h5ad_input)
        annotation_dict = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', annotation_dict)

# %%
adata = sc.read_h5ad(h5ad_input)
mapping_df = pd.read_csv(annotation_dict)

# %%
celltype_mapping = dict(zip(mapping_df['celltype1'],
                                       mapping_df['LR_Cluster']))
print("INFO: mapping")
print(celltype_mapping)

celltype_mapping2 = dict(zip(mapping_df['celltype1'],
                                       mapping_df['LR_Cluster2']))
print(celltype_mapping)

# %%
adata.obs['LR_Cluster'] = adata.obs['celltype'].map(celltype_mapping)
adata.obs['LR_Cluster2'] = adata.obs['celltype'].map(celltype_mapping2)

print("INFO: adata.obs")
print(adata.obs.head())

# %%
adata.write_h5ad(h5ad_output)