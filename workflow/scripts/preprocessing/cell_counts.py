# (C) Robin Fallegger
# %%
import scanpy as sc
import pandas as pd
import platform as platform
import os as os

# %%
# if snakemake object is in local
if 'snakemake' in locals():
    input_file = snakemake.input[0]
else:
    wildcards = {'dataset': 'Julio18'}
    input_file = 'sd22b002/projects/MichiganIgAN/results/preprocessing/{dataset}.h5ad'.format(dataset=wildcards['dataset'])

    if platform.system() == 'Linux':
        input_file = os.path.join('/mnt/sds-hd/', input_file)
    elif platform.system() == 'Darwin':
        input_file = os.path.join('/Volumes', input_file)

# %%
#check if input file exists
if not os.path.exists(input_file):
    raise FileNotFoundError(input_file + ' does not exist')

### Load RNA modalities
print('INFO: Loading assays')
adata = sc.read_h5ad(input_file)
print(adata.obs)
print(list(adata.obs))
# %%
# keep only relevant columns
print('INFO: Filtering adata.obs object')
cols_to_keep = ['ID', 'celltype']
cols_to_keep = adata.obs.columns.intersection(cols_to_keep)
cells = adata.obs.filter(cols_to_keep, axis='columns')

# check that specific columns were kept
# TODO: make a function to choose these automatically based on available columns
# e.g. hybrid cell annotation if there and else celltype_2
obligatory_cols = ['ID', 'celltype']
assert all([col in cells.columns for col in obligatory_cols ]), 'Not all columns were kept'

# filter out NA cells
cells = cells[cells['celltype'] != 'NA']
print(cells)
print(cells.shape)
print(cols_to_keep)
print([*cols_to_keep])

#print('INFO: Counting cells')
# %% get cells counts per celltype per patient, and the cellstate proportions
cell_count = cells.groupby([*cols_to_keep], observed=True).size().rename('n_cells').reset_index()
cell_count = cell_count.assign(total_cells = lambda x: x.groupby(['ID', 'celltype'])['n_cells'].transform(lambda x: x.sum()),
            cs_prop = lambda x: x['n_cells']/x['total_cells'])

# %%

if 'snakemake' in locals():
    output_file = snakemake.output[0]
    cell_count.to_csv(output_file, sep=',', index=False)