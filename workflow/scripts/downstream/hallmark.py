# 15.08.2024 Charlotte Boys
# %% Import libraries
import platform as platform
import os as os
import decoupler as dc
import pandas as pd
import numpy as np

# %% Parse snakemake object
if 'snakemake' in locals():
    loadings_path = snakemake.input['loadings_path']
    acts = snakemake.output['acts']
else:
    import os
    os.chdir('/Users/charlotteboys/GitHub/OCEAN_MCFA')
    loadings_path = 'results/mofacell/loadings~combined_sn10xLRCluster2.csv'


# %%
loadings_all = pd.read_csv(loadings_path, index_col=0)
loadings = loadings_all[~loadings_all['view'].str.contains('&')]
cell_types = np.unique(loadings['view'])

# %% 
msigdb = dc.get_resource('MSigDB')

hallmark = msigdb[msigdb['collection']=='hallmark']
# Remove duplicated entries
hallmark = hallmark[~hallmark.duplicated(['geneset', 'genesymbol'])]

# %%
def enrichment(factor,cell_type, net, positive = True, ascending=False, n=200):
    # Ascending = False -> positive loadings
    # Ascending = True -> negative loadings
    loadings_df = loadings[
    (loadings['Factor'] == factor) &
    (loadings['view'] == cell_type)
].sort_values('Loading', ascending=not(positive)).head(n)
    enr = dc.get_ora_df(
        df = loadings_df['variable'].to_numpy(),
        net = net,
        source = 'geneset',
        target = 'genesymbol'
    )
    return enr.sort_values('FDR p-value')

# %%   
enrichment(factor='Factor3',
           cell_type = 'DTL_aPT',
            net = hallmark,
            positive = True).head(15)


# %%
positive_or_not =  [True, False]

# Initialize lists to store results
enr_list = []

# Run pathway_activity for each factor and cell type
for direction in positive_or_not:
    results = [enrichment(factor = 'Factor3',
                          cell_type=cell_type,
                            net = hallmark,
                            positive = direction).assign(view = cell_type) for cell_type in cell_types]
    
    # Add direction column to the results
    enr_with_dir = pd.concat(
        [result.assign(positive=direction) for result in results],
                                                 axis=0)
    
    enr_list.append(enr_with_dir)

# Combine all results into a single DataFrame
enr_combined = pd.concat(enr_list, axis=0)

enr_combined.to_csv('results/downstream/hallmark_factor3.csv')


# %%
# Initialize lists to store results
enr_list = []

# Run pathway_activity for each factor and cell type
for direction in positive_or_not:
    results = [enrichment(factor = 'Factor5',
                          cell_type=cell_type,
                            net = hallmark,
                            positive = direction).assign(view = cell_type) for cell_type in cell_types]
    
    # Add direction column to the results
    enr_with_dir = pd.concat(
        [result.assign(positive=direction) for result in results],
                                                 axis=0)
    
    enr_list.append(enr_with_dir)

# Combine all results into a single DataFrame
enr_combined = pd.concat(enr_list, axis=0)

enr_combined.to_csv('results/downstream/hallmark_factor5.csv')