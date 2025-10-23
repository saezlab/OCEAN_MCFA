rule create_mdata: # Done
    input:
        h5ad_input = 'results/preprocessing/{dataset}.h5ad',
        coldata_path = 'data/metadata_clean.csv'#'data/Metadata_For_Charlotte.csv'
    output:
        h5ad_output = 'results/mofacell/{dataset}~{subset}.h5ad',
        h5mu_output = 'results/mofacell/{dataset}~{subset}.h5mu'
    resources:
        mem_mb=100000
    wildcard_constraints:
        subset = 'sn10xcanon2'
    params:
        sample_key = 'EdgarID',
        condition_key = 'Cohort',
        min_cells = 50,
        batch_key = 'Project1',
        groupby = lambda wildcards: config["subsets"][wildcards.subset]["groupby"],
        sample_subset = lambda wildcards: config["subsets"][wildcards.subset]["sample_subset"],
        cell_subset = lambda wildcards: config["subsets"][wildcards.subset]["cell_subset"],
        drop_views = lambda wildcards: config["mofacell"][wildcards.dataset][wildcards.subset]["drop_views"]
    singularity:
        'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/mofacell/create_mdata.py'

rule run_liana:
    input:
        h5ad_input = 'results/mofacell/{dataset}~{subset}.h5ad'
    output:
        h5ad_output = 'results/mofacell/LIANA~{dataset}~{subset}.h5ad'
    resources:
        mem_mb=100000
    params:
        sample_key = 'EdgarID',
        groupby = lambda wildcards: config["subsets"][wildcards.subset]["groupby"]
    singularity:
        'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/mofacell/run_liana.py'

rule create_liana_mdata:
    input:
        h5ad_input = 'results/mofacell/LIANA~{dataset}~{subset}.h5ad'
    output:
        h5mu_output = 'results/mofacell/LIANA~{dataset}~{subset}.h5mu'
    resources:
        mem_mb=100000
    params:
        sample_key = 'EdgarID',
        groupby = lambda wildcards: config["subsets"][wildcards.subset]["groupby"]
    singularity:
        'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/mofacell/create_liana_mdata.py'

rule combine_liana_mdata:
    input:
        h5mu_liana = lambda wildcards: config["niches"][wildcards.dataset][wildcards.niche]["h5mu_liana"],
        h5mu_input = lambda wildcards: config["niches"][wildcards.dataset][wildcards.niche]["h5mu_input"]
    output:
        h5mu_output = 'results/mofacell/{dataset}~LIANA_{niche}.h5mu'
    wildcard_constraints:
        niche = 'all'
    params:
        cell_subset = lambda wildcards: config["niches"][wildcards.dataset][wildcards.niche]["cell_subset"],
        hvg = 'hvg'
    singularity:
        'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/mofacell/combine_liana_mdata.py'

rule run_mofa_subset: # Done
    input:
        h5mu_input = 'results/mofacell/{dataset}~{subset}.h5mu'
    output:
        model_output = 'results/mofacell/model~{dataset}~{subset}~{hvg}.h5ad',
        scores_output = 'results/mofacell/scores~{dataset}~{subset}~{hvg}.csv',
        scores_wide_output = 'results/mofacell/widescores~{dataset}~{subset}~{hvg}.csv',
        loadings_output = 'results/mofacell/loadings~{dataset}~{subset}~{hvg}.csv'
    wildcard_constraints:
        hvg = 'all|hvg'
    params:
        nfactors = 20,
        hvg = lambda wildcards: wildcards.hvg,
        sample_key = 'EdgarID'
    singularity:
        'workflow/envs/mofaxmetacells.0.0.2.sif'
    script:
        '../scripts/mofacell/run_mofa_subset.py'

rule calculate_r2_subset: # Done, updated 21.02.24
    input:
        model_input = 'results/mofacell/model~{dataset}~{subset}~{hvg}.h5ad'
    output:
        variance_explained_by_sample_output = 'results/mofacell/r2_by_sample~{dataset}~{subset}~{hvg}.csv',
        variance_explained_total_output = 'results/mofacell/r2_total~{dataset}~{subset}~{hvg}.csv'
    params:
        condition_key = 'Cohort'
    singularity:
        'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/mofacell/calculate_r2.py'

rule plot_r2_scores: # Done 21.02.24
    input:
        r2 = 'results/mofacell/r2_by_sample~{dataset}~{subset}~{hvg}.csv',
        scores = 'results/mofacell/scores~{dataset}~{subset}~{hvg}.csv'
    output:
        plt = 'plots/mofacell/scores_r2_boxplot~{dataset}~{subset}~{hvg}.pdf',
        ggs = 'results/mofacell/scores_r2_boxplot~{dataset}~{subset}~{hvg}.rds'
    singularity:
        'workflow/envs/R_env.sif'
    script:
        '../scripts/mofacell/plot_r2_scores.R'

rule plot_r2_scores_by_project: # Done 12.03.24
    input:
        r2 = 'results/mofacell/r2_by_sample~{dataset}~{subset}~{hvg}.csv',
        scores = 'results/mofacell/scores~{dataset}~{subset}~{hvg}.csv',
        metadata = 'data/Metadata_For_Charlotte.csv',
        r2_total = 'results/mofacell/r2_total~{dataset}~{subset}~{hvg}.csv'
    output:
        plt = 'plots/mofacell/scores_r2_boxplot_by_project~{dataset}~{subset}~{hvg}.pdf',
        ggs = 'results/mofacell/scores_r2_boxplot_by_project~{dataset}~{subset}~{hvg}.rds'
    singularity:
        'workflow/envs/R_env.sif'
    script:
        '../scripts/mofacell/plot_r2_scores_by_project.R'

rule plot_r2_scores_by_disease: # Done 12.03.24
    input:
        r2 = 'results/mofacell/r2_by_sample~{dataset}~{subset}~{hvg}.csv',
        scores = 'results/mofacell/scores~{dataset}~{subset}~{hvg}.csv',
        r2_total = 'results/mofacell/r2_total~{dataset}~{subset}~{hvg}.csv'
    output:
        plt = 'plots/mofacell/scores_r2_boxplot_by_disease~{dataset}~{subset}~{hvg}.pdf',
        ggs = 'results/mofacell/scores_r2_boxplot_by_disease~{dataset}~{subset}~{hvg}.rds'
    singularity:
        'workflow/envs/R_env.sif'
    script:
        '../scripts/mofacell/plot_r2_scores_by_disease.R'
    

rule plot_r2_total:
    input:
        r2 = 'results/mofacell/r2_total~{dataset}~{subset}~{hvg}.csv'
    output:
        data = 'results/mofacell/r2_total_heatmap~{dataset}~{subset}~{hvg}.csv',
        plt = 'plots/mofacell/r2_total_heatmap~{dataset}~{subset}~{hvg}.pdf',
        gg = 'results/mofacell/r2_total_heatmap~{dataset}~{subset}~{hvg}.rds'
    singularity:
        'workflow/envs/R_env.sif'
    script:
        '../scripts/mofacell/plot_r2_total.R'

rule test_associations:
    input:
        scores_wide = 'results/mofacell/widescores~{dataset}~{subset}~{hvg}.csv',
        metadata_csv = 'results/preprocessing/metadata_clean.csv',
        clr = lambda wildcards: config["clr"][wildcards.subset],
        r2_csv = 'results/mofacell/r2_by_sample~{dataset}~{subset}~{hvg}.csv',
        r2_total_csv = 'results/mofacell/r2_total~{dataset}~{subset}~{hvg}.csv'
    output:
        associations_csv = 'results/mofacell/associations~{dataset}~{subset}~{hvg}.csv',
        factors_obs_csv = 'results/mofacell/factors_obs~{dataset}~{subset}~{hvg}.csv'
    singularity:
        'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/mofacell/test_associations.py'

rule plot_associations:
    input:
        associations = 'results/mofacell/associations~{dataset}~{subset}~{hvg}.csv',
    output:
        plt1 = 'plots/mofacell/associations_clr~{dataset}~{subset}~{hvg}.pdf',
        plt1png = 'plots/mofacell/associations_clr~{dataset}~{subset}~{hvg}.png',
        plt2 = 'plots/mofacell/associations_cov~{dataset}~{subset}~{hvg}.pdf',
        plt2png = 'plots/mofacell/associations_cov~{dataset}~{subset}~{hvg}.png',
    singularity:
        'workflow/envs/R_env.sif'
    script:
        '../scripts/mofacell/plot_associations.R'
#######################
#
# QC of data for MOFAcell models
#
#######################
# Here are rules to plot the plot some QC metrics on the data used to plot MOFAcell models
# in particular:
# - plot the number of cells, genes and counts per sample

rule plot_qc:
    input:
        'results/mofacell/{dataset}~{subset}.h5mu'
        #'results/mofacell/{omic}/{source}/{psbulk}/{proj_omic}_projected.h5mu'
    params:
        colors = config['general_plotting']['colors'],
        sample_key = 'ID',
        condition_key = 'Cohort',
        omic = lambda wildcards: config["subsets"][wildcards.subset]["omic"]
    output:
        'plots/mofacell/{dataset}~{subset}_overall_qc.pdf' #{omic}/{source}/{psbulk}/overall_qc_w_{proj_omic}.pdf'
    singularity:
        'workflow/envs/scanpy.0.0.3.sif'
    script:
        '../scripts/mofacell/plot_mofa_qc.py'