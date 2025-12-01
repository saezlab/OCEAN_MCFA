rule clr:
    input:
        h5ad_input = config["input"]
    output:
        clr_csv = 'results/preprocessing/clr.csv',
        prop_csv = 'results/preprocessing/prop.csv'
    params:
        groupby = config["clr"]["groupby"], #'LR_Cluster',
        sample_key = config['sample_key'],
    resources:
        mem_mb=50000
    singularity:
        config["singularities"]["scanpy"]#'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/preprocessing/clr.py'
