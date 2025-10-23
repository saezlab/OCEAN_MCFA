rule run_mofa_robustness: # Done
    input:
        h5mu_input = ancient('results/mofacell/{dataset}~{subset}.h5mu')
    output:
        model_output = 'results/robustness/model~{dataset}~{subset}~{hvg}~{nfactors}.h5ad',
        scores_output = 'results/robustness/scores~{dataset}~{subset}~{hvg}~{nfactors}.csv',
        scores_wide_output = 'results/robustness/widescores~{dataset}~{subset}~{hvg}~{nfactors}.csv',
        loadings_output = 'results/robustness/loadings~{dataset}~{subset}~{hvg}~{nfactors}.csv'
    params:
        nfactors = lambda wildcards: wildcards.nfactors 
        hvg = lambda wildcards: wildcards.hvg,
        sample_key = 'EdgarID'
    wildcard_constraints:
        hvg = 'all|hvg'
    singularity:
        'workflow/envs/mofaxmetacells.0.0.2.sif'
    script:
        '../scripts/mofacell/run_mofa_subset.py'

rule calculate_r2_robustness: # Done, updated 21.02.24
    input:
        model_input = 'results/robustness/model~{dataset}~{subset}~{hvg}~{nfactors}.h5ad'
    output:
        variance_explained_by_sample_output = 'results/robustness/r2_by_sample~{dataset}~{subset}~{hvg}~{nfactors}.csv',
        variance_explained_total_output = 'results/robustness/r2_total~{dataset}~{subset}~{hvg}~{nfactors}.csv'
    params:
        condition_key = 'Cohort'
    singularity:
        'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/mofacell/calculate_r2.py'

rule combine_scores:
    input:
        expand('results/robustness/scores~{{dataset}}~{{subset}}~{{hvg}}~{nfactors}.csv', nfactors=config["nfactors_list"])
    output:
        'results/robustness/scores~{dataset}~{subset}~{hvg}.csv'
    singularity:
        'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/robustness/combine.py'

rule combine_loadings:
    input:
        expand('results/robustness/loadings~{{dataset}}~{{subset}}~{{hvg}}~{nfactors}.csv', nfactors=config["nfactors_list"])
    output:
        'results/robustness/loadings~{dataset}~{subset}~{hvg}.csv'
    singularity:
        'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/robustness/combine.py'

rule combine_r2_by_sample:
    input:
        expand('results/robustness/r2_by_sample~{{dataset}}~{{subset}}~{{hvg}}~{nfactors}.csv', nfactors=config["nfactors_list"])
    output:
        'results/robustness/r2_by_sample~{dataset}~{subset}~{hvg}.csv'
    singularity:
        'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/robustness/combine.py'

rule combine_r2_total:
    input:
        expand('results/robustness/r2_total~{{dataset}}~{{subset}}~{{hvg}}~{nfactors}.csv', nfactors=config["nfactors_list"])
    output:
        'results/robustness/r2_total~{dataset}~{subset}~{hvg}.csv'
    singularity:
        'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/robustness/combine.py'