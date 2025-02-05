# IMPORTS

import pandas

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DEFINE VARIABLES

# Define the data directory, explicitly
# NOTE - this is where the outputs of the pipeline processing will go into for later analysis
data_dir = '/data/CARD_singlecell/Brain_atlas/SN_Multiome/'

# Define the working directory, explictly
work_dir = '/data/CARD_singlecell/SN_control_atlas/multiome-pipeline/'

# Number of threads to use when running the rules
num_workers = 8

# Define where the metadata data exists for each sample to be processed
input_table = '/data/CARD_singlecell/SN_control_atlas/input/SN_control_samples.csv'

# Read in the list of 
batches = pandas.read_csv(input_table)['Use_batch'].tolist()
samples = pandas.read_csv(input_table)['Sample'].tolist()

# Define environmental .YML files as a dictionary
envs = {
    'cellbender'        : 'envs/cellbender.yml', 
    # 'muon'              : 'envs/muon.yml',
    'singlecell'        : 'envs/single_cell_cpu.yml', 
    'single_cell_gpu'   : 'envs/single_cell_gpu.yml', 
    'atac'              : 'envs/snapATAC2.yml'
    }

# Define RNA thresholds
mito_percent_thresh = 15
ribo_percent_thresh = 10
doublet_thresh      = 0.15
min_genes_per_cell  = 250

# Define ATAC thresholds
min_peak_counts         = 500
min_num_cell_by_counts  = 10

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SNAKEMAKE RULES

rule all:
    input:
        rna_anndata = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad', 
            zip,
            batch = batches,
            sample = samples
            )

"""
# NOTE - only needs to run once at the beginning

rule cellbender:
    script:
        work_dir+'scripts/cellbender_array.sh'

rule pileup:
    script:
        'scripts/pileup.sh'
"""   

# ADDING GVCF, QTL, work here

rule preprocess:
    input:
        input_table = input_table,
        rna_anndata = data_dir+'batch{batch}/Multiome/{dataset}-ARC/outs/cellbender_gex_counts_filtered.h5'
    output:
        rna_anndata = data_dir+'batch{batch}/Multiome/{dataset}-ARC/outs/01_{dataset}_anndata_object_rna.h5ad'
    conda:
        envs['singlecell']
    params:
        sample = '{dataset}'
    resources:
        runtime=120, mem_mb=64000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'scripts/preprocess.py'

rule merge_unfiltered:
    input:
        rna_anndata = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad', 
            zip,
            batch = batches,
            sample = samples
            )
    output:
        merged_rna_anndata = data_dir+'control_atlas/01_merged_anndata_rna.h5ad'
    conda:
        envs['singlecell']
    params:
        samples = samples
    resources:
        runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'scripts/merge_anndata.py'

rule plot_qc_rna:
    input:
        merged_rna_anndata = data_dir+'control_atlas/01_merged_anndata_rna.h5ad'
    conda:
        # envs['muon']
        envs['singlecell']
    resources:
        runtime=960, mem_mb=500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'scripts/plot_qc_metrics.py'

rule filter_rna:
    input:        
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad'
    output:
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/02_{sample}_anndata_filtered_rna.h5ad'
    conda:
        # envs['muon']
        envs['singlecell']
    params:
        mito_percent_thresh     = mito_percent_thresh,
        doublet_thresh          = doublet_thresh,
        min_genes_per_cell      = min_genes_per_cell,
        ribo_percent_thresh     = ribo_percent_thresh
    resources:
        runtime=120, mem_mb=100000, disk_mb=10000, slurm_partition='quick' 
    script: 
        work_dir+'scripts/rna_filter.py'

rule merge_filtered_rna:
    input:
        rna_anndata = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/02_{sample}_anndata_filtered_rna.h5ad', 
            zip,
            batch = batches,
            sample = samples
            )
    output:
        merged_rna_anndata = data_dir+'control_atlas/02_filtered_anndata_rna.h5ad'
    conda:
        envs['singlecell']
    params:
        samples = samples
    resources:
        runtime=120, mem_mb=1000000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'scripts/merge_anndata.py'

rule atac_preprocess:
    input:
        fragment_file=data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz'
    output:
        atac_anndata=data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_atac-control.h5ad'
    conda:
        envs['atac']
    resources:
        runtime=120, mem_mb=50000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'scripts/atac_preprocess.py'

"""
# Get back to this once separate script written
rule merge_unfiltered_atac:
    input:
        atac_anndata=expand(
            work_dir+'data/samples/{dataset}/01_{dataset}_anndata_object_atac-control.h5ad', 
            dataset=datasets
            )
    output:
        merged_atac_anndata=work_dir+'data/control_atlas/01_merged_anndata_atac-control.h5ad'
    conda:
        envs['singlecell']
    resources:
        runtime=480, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'scripts/merge_atac.py'
"""

rule plot_qc_atac:
    input:
        atac_anndata=data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad'
    conda:
        envs['atac']
    resources:
        runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem'
    script:
        work_dir+'scripts/atac_plot_qc.py'

rule filter_atac:
    input:
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/02_{sample}_anndata_filtered_rna.h5ad',
        atac_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_atac-control.h5ad'
    output:
        atac_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac-control.h5ad',
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_filtered_rna.h5ad'
    conda:
        envs['atac']
    resources:
        runtime=30, mem_mb=50000, slurm_partition='quick'
    script:
        work_dir+'scripts/atac_filter.py'

rule merge_multiome_rna:
    input:
        rna_anndata = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_filtered_rna.h5ad', 
            zip,
            batch = batches,
            sample = samples
            )
    output:
        merged_rna_anndata = data_dir+'control_atlas/03_filtered_anndata_rna.h5ad'
    conda:
        envs['singlecell']
    params:
        samples = samples
    resources:
        runtime=120, mem_mb=1000000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'scripts/merge_anndata.py'

"""
rule rna_atac_filter:
    input:
        merged_rna_anndata = work_dir+'data/control_atlas/02_filtered_anndata_rna.h5ad',
        merged_atac_anndata = work_dir+'data/control_atlas/02_filtered_anndata_atac-control.h5ad'
    output:
        merged_rna_anndata = work_dir+'data/control_atlas/03_filtered_anndata_rna.h5ad',
        merged_atac_anndata = work_dir+'data/control_atlas/03_filtered_anndata_atac-control.h5ad'
    conda:
        envs['singlecell']
    resources:
        runtime=2880, mem_mb=3000000, slurm_partition='largemem'
    script:
        'scripts/filter_rna_atac.py'
"""

rule rna_model:
    input:
        merged_rna_anndata = data_dir+'control_atlas/03_filtered_anndata_rna.h5ad'
    output:
        merged_rna_anndata = data_dir+'control_atlas/04_modeled_anndata_rna.h5ad',
        model_history = work_dir+'model_elbo/rna_model_history.csv'
    params:
        model = work_dir+'data/models/rna/'
    conda:
        envs['single_cell_gpu']
    threads:
        64
    resources:  
        runtime=2880, disk_mb=500000, mem_mb=300000, gpu=2, gpu_model='v100x'
    script:
        'scripts/rna_model.py'

rule annotate:
    input:
        merged_rna_anndata = data_dir+'control_atlas/04_modeled_anndata_rna.h5ad'
    output:
        merged_rna_anndata = data_dir+'control_atlas/05_annotated_anndata_rna.h5ad',
        cell_annotate = work_dir+'data/rna_cell_annot.csv'
    conda:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=500000, slurm_partition='largemem'
    script:
        'scripts/annotate.py'

rule SCANVI_annot:
    input:
        merged_rna_anndata = data_dir+'control_atlas/05_annotated_anndata_rna.h5ad',
        model = work_dir+'data/models/rna'
    output:
        merged_rna_anndata = data_dir+'control_atlas/06_SCVI_anndata_rna.h5ad',
        model_history = work_dir+'model_elbo/SCANVI_model_history.csv'
    conda:
        envs['single_cell_gpu']
    threads:
        64
    resources:
        runtime=2880, disk_mb=500000, mem_mb=300000, gpu=2, gpu_model='v100x'
    script:
        'scripts/SCANVI_annot.py'

rule merge_multiome_atac:
    input:
        cell_annotate = data_dir+'data/rna_cell_annot.csv',
        atac_anndata = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac-control.h5ad', 
            zip,
            batch = batches,
            sample = samples
            )
    output:
        umap_data = data_dir+'data/atac_umap.csv',
        var_data = data_dir+'data/atac_var_selected.csv',
        merged_atac_anndata = data_dir+'control_atlas/03_filtered_anndata_atac-control.h5ad'
    conda:
        envs['atac']
    resources:
        runtime=2880, mem_mb=3000000, slurm_partition='largemem'
    threads:
        64
    script:
        work_dir+'scripts/merge_atac.py' 

"""
rule atac_model:
    input:
        atac_anndata = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac-control.h5ad', 
            batch = batches,
            sample = samples
            )
    output:
        atac_umap = work_dir+'data/atac_umap.csv',
        atac_var = work_dir+'data/atac_var_selected.csv'
    conda:
        envs['atac']
    threads:
        64
    resources:
        runtime=2880, disk_mb=500000, mem_mb=300000, gpu=4
    script:
        'scripts/atac_model.py'

rule apply_atac_annotation:
rule DGE:
rule DA:
"""