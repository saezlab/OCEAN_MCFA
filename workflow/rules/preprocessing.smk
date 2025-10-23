# # Complete 15.02.24
# rule rds_to_seurat:
#     input:
#         data = 'data/{dataset}.rds'
#     output:
#         seurat = 'results/preprocessing/temp/{dataset}.h5Seurat'
#     resources:
#         mem_mb=100000
#     singularity:
#         'workflow/envs/dot.0.0.1.sif'
#     script:
#         "../scripts/preprocessing/rds_to_h5Seurat.R"
# # Complete 15.02.24
# rule seurat_to_h5ad:
#     input:
#         data = 'results/preprocessing/temp/{dataset}.h5Seurat'
#     output:
#         ad = 'results/preprocessing/temp/{dataset}.h5ad'
#     resources:
#         mem_mb=100000
#     wildcard_constraints:
#         dataset='Julio_OCEAN_Nereid'
#     params:
#         assays = lambda w: config['preprocessing'][w.dataset].get('assays', {}),
#         reductions = lambda w: config['preprocessing'][w.dataset].get('reductions', {})
#     singularity:
#         'workflow/envs/scanpy.0.0.5.sif'
#     script:
#         "../scripts/preprocessing/h5Seurat_to_h5ad.py"

rule add_annotations:
    input:
        h5ad_input = 'data/{dataset}.h5ad'#'results/preprocessing/temp/{dataset}.h5ad',
        annotation_dict = 'data/OCEAN_LR.csv'
    output:
        h5ad_output = 'results/preprocessing/{dataset}.h5ad'
    resources:
        mem_mb=100000
    wildcard_constraints:
        dataset='Julio_OCEAN_Nereid'
    singularity:
        'workflow/envs/scanpy.0.0.5.sif'
    script:
        "../scripts/preprocessing/add_annotations.py"

rule clr:
    input:
        h5ad_input = 'results/preprocessing/{dataset}.h5ad'
    output:
        clr_csv = 'results/preprocessing/clr~{dataset}~{subset}.csv',
        prop_csv = 'results/preprocessing/prop~{dataset}~{subset}.csv'
    # wildcard_constraints:
    #     subset = 'all|sn10xcanon|sn10xcanon2|sn10x|sn|glomerular'
    params:
        groupby = lambda wildcards: config["subsets"][wildcards.subset]["groupby"],
        cell_subset = 'all', #lambda wildcards: config["subsets"][wildcards.subset]["cell_subset"],
        sample_key = 'EdgarID'
    resources:
        mem_mb=50000
    singularity:
        'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/preprocessing/clr.py'

# rule clean_up_metadata:
#     input:
#         metadata = 'data/Metadata_For_Charlotte.csv'
#     output:
#         metadata_clean = 'results/preprocessing/metadata_clean.csv'
#     singularity:
#         'workflow/envs/R_env.sif'
#     script:
#         '../scripts/preprocessing/clean_up_metadata.R'
