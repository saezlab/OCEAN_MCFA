rule create_mdata:
    input:
        h5ad_input = config["input"],
        coldata_path = config["metadata"],
        annotation_dict = 'annotations/OCEAN_LR.csv'
    output:
        h5ad_output = 'results/mofacell/{model}.h5ad',
        h5mu_output = 'results/mofacell/{model}.h5mu'
    resources:
        mem_mb=100000
    wildcard_constraints:
        model = 'sn10xLRCluster|sn10xLRCluster2'
    params:
        sample_key = config['sample_key'],
        batch_key = config['batch_key'],
        condition_key = config['condition_key'],
        groupby = lambda wildcards: config["model"][wildcards.model]["groupby"],
        min_cells = lambda wildcards: config["model"][wildcards.model]['min_cells'],
        sample_subset = lambda wildcards: config["model"][wildcards.model]["sample_subset"],
        drop_views = lambda wildcards: config["model"][wildcards.model]["drop_views"]
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/mofacell/create_mdata.py'

rule run_liana:
    input:
        h5ad_input = 'results/mofacell/{model}.h5ad',
    output:
        h5ad_output = 'results/mofacell/LIANA_{model}.h5ad'
    resources:
        mem_mb=100000
    params:
        sample_key = config['sample_key'],
        groupby = lambda wildcards: config["model"][wildcards.model]["groupby"]
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/mofacell/run_liana.py'

rule create_liana_mdata:
    input:
        h5ad_input = 'results/mofacell/LIANA_{model}.h5ad'
    output:
        h5mu_output = 'results/mofacell/LIANA_{model}.h5mu'
    resources:
        mem_mb=100000
    params:
        sample_key = config['sample_key'],
        batch_key = config['batch_key'],
        condition_key = config['condition_key'],
        groupby = lambda wildcards: config["model"][wildcards.model]["groupby"]
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/mofacell/create_liana_mdata.py'

rule combine_liana_mdata:
    input:
        h5mu_liana = 'results/mofacell/LIANA_{model}.h5mu',
        h5mu_input = 'results/mofacell/{model}.h5mu'
    output:
        h5mu_output = 'results/mofacell/combined_{model}.h5mu'
    params:
       hvg = lambda wildcards: config["model"][wildcards.model]['hvg']
    wildcard_constraints:
        hvg = 'all|hvg'
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/mofacell/combine_liana_mdata.py'

rule run_mofa_subset:
    input:
        h5mu_input = 'results/mofacell/combined_{model}.h5mu'
    output:
        model_output = 'results/mofacell/model~combined_{model}.h5ad',
        scores_output = 'results/mofacell/scores~combined_{model}.csv',
        scores_wide_output = 'results/mofacell/widescores~combined_{model}.csv',
        loadings_output = 'results/mofacell/loadings~combined_{model}.csv'
    params:
        nfactors = lambda wildcards: config["model"][wildcards.model]['nfactors'],
        hvg = lambda wildcards: config["model"][wildcards.model]['hvg'],
        sample_key = config['sample_key'],
    wildcard_constraints:
        hvg = 'all|hvg',
    singularity:
        config["singularities"]["mofaxmetacells"]
    script:
        '../scripts/mofacell/run_mofa_subset.py'

rule calculate_r2_subset: # Done, updated 21.02.24
    input:
        model_input = 'results/mofacell/model~combined_{model}.h5ad'
    output:
        variance_explained_by_sample_output = 'results/mofacell/r2_by_sample~combined_{model}.csv',
        variance_explained_total_output = 'results/mofacell/r2_total~combined_{model}.csv'
    params:
        condition_key = config['condition_key'],
        sample_key = config['sample_key']
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/mofacell/calculate_r2.py'