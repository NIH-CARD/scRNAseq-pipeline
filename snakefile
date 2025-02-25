import pandas as pd
import os

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""

# Define the data directory, explicitly
data_dir = '/data/CARD_singlecell/Brain_atlas/SN_Multiome/'
# Define the working directory, explictly as the directory of this pipeline
work_dir = os.getcwd()

# Number of threads to use when running the rules
num_workers = 8

# Define where the metadata data exists for each sample to be processed
metadata_table = work_dir+'/input/example_metadata.csv'
# Define where celltypes/cell marker gene 
gene_markers_file = work_dir+'/input/example_marker_genes.csv'

# Read in the list of batches and samples
batches = pd.read_csv(metadata_table)['Use_batch'].tolist()
samples = pd.read_csv(metadata_table)['Sample'].tolist()

# Name of the disease parameter
disease_param = 'Primary Diagnosis'
# Define disease states
control = 'control'
diseases = ['PD', 'DLB']

# Define the cell types to look for, from gene marker file
cell_types = pd.read_csv(gene_markers_file)['cell type']

# Define RNA thresholds
mito_percent_thresh = 15
ribo_percent_thresh = 10
doublet_thresh = 0.15
min_genes_per_cell = 250

# Define ATAC thresholds
min_peak_counts = 500
min_num_cell_by_counts = 10

"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Singularity containers to be downloaded from Quay.io, done in snakemake.sh
envs = {
    'singlecell': 'envs/single_cell_cpu.sif', 
    'atac': 'envs/snapATAC2.sif'
    }

rule all:
    input:
        rna_anndata=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac.h5ad', 
            zip,
            batch=batches,
            sample=samples
            ),
        atac_anndata = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac.h5ad',
            zip,
            sample=samples,
            batch=batches
            ),
        merged_multiome = work_dir+'/atlas/multiome_atlas.h5mu',
        output_DGE_data = expand(
            work_dir + '/data/significant_genes/rna/rna_{cell_type}_{disease}_DGE.csv',
            cell_type = cell_types,
            disease = diseases
            ),
        output_DAR_data = expand(
            work_dir + '/data/significant_genes/atac/atac_{cell_type}_{disease}_DAR.csv',
            cell_type = cell_types,
            disease = diseases
            )
        
# This needs to be forced to run once
rule cellbender:
    script:
        work_dir+'/scripts/cellbender_array.sh'

rule preprocess:
    input:
        metadata_table=metadata_table,
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/cellbender_gex_counts_filtered.h5'
    output:
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        sample='{sample}'
    resources:
        runtime=120, mem_mb=64000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'/scripts/preprocess.py'

rule merge_unfiltered:
    input:
        rna_anndata=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        merged_rna_anndata = work_dir+'/atlas/01_merged_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        samples=samples
    resources:
        runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_anndata.py'

rule plot_qc_rna:
    input:
        merged_rna_anndata = work_dir+'/atlas/01_merged_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    resources:
        runtime=960, mem_mb=500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/plot_qc_metrics.py'

rule filter_rna:
    input:        
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad'
    output:
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/02_{sample}_anndata_filtered_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        mito_percent_thresh = mito_percent_thresh,
        doublet_thresh = doublet_thresh,
        min_genes_per_cell = min_genes_per_cell,
        ribo_percent_thresh = ribo_percent_thresh
    resources:
        runtime=120, mem_mb=100000, disk_mb=10000, slurm_partition='quick' 
    script: 
        work_dir+'/scripts/rna_filter.py'

rule merge_filtered_rna:
    input:
        rna_anndata=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/02_{sample}_anndata_filtered_rna.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        merged_rna_anndata = work_dir+'/atlas/02_filtered_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        samples=samples
    resources:
        runtime=120, mem_mb=1000000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_anndata.py'

rule atac_preprocess:
    input:
        fragment_file=data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz'
    output:
        atac_anndata=data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_atac.h5ad'
    singularity:
        envs['atac']
    resources:
        runtime=120, mem_mb=50000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'/scripts/atac_preprocess.py'

rule merge_unfiltered_atac:
    input:
        rna_anndata=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_atac.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        merged_atac_anndata = work_dir+'/atlas/01_merged_anndata_atac.h5ad'
    singularity:
        envs['singlecell']
    resources:
        runtime=480, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_atac.py'

rule plot_qc_atac:
    input:
        atac_anndata=data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad'
    singularity:
        envs['atac']
    resources:
        runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem'
    script:
        work_dir+'/scripts/atac_plot_qc.py'

rule filter_atac:
    input:
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/02_{sample}_anndata_filtered_rna.h5ad',
        atac_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_atac.h5ad'
    output:
        atac_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac.h5ad',
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_filtered_rna.h5ad'
    singularity:
        envs['atac']
    resources:
        runtime=30, mem_mb=50000, slurm_partition='quick'
    script:
        work_dir+'/scripts/atac_filter.py'

rule merge_multiome_rna:
    input:
        rna_anndata=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_filtered_rna.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        merged_rna_anndata = work_dir+'/atlas/03_filtered_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        samples=samples
    resources:
        runtime=120, mem_mb=1000000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_anndata.py'

rule rna_model:
    input:
        merged_rna_anndata = work_dir+'/atlas/03_filtered_anndata_rna.h5ad'
    output:
        merged_rna_anndata = work_dir+'/atlas/04_modeled_anndata_rna.h5ad',
        model_history = work_dir+'/model_elbo/rna_model_history.csv'
    params:
        model = work_dir+'/data/models/rna/'
    singularity:
        envs['singlecell'] # GPU environment needs work: envs['single_cell_gpu']
    threads:
        64
    resources:
        runtime=2880, disk_mb=500000, mem_mb=300000#, gpu=2, gpu_model='v100x'
    script:
        'scripts/rna_model.py'

rule annotate:
    input:
        merged_rna_anndata = work_dir+'/atlas/04_modeled_anndata_rna.h5ad',
        gene_markers = gene_markers_file
    output:
        merged_rna_anndata = work_dir+'/atlas/05_annotated_anndata_rna.h5ad',
        cell_annotate = work_dir+'/data/rna_cell_annot.csv'
    singularity:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=500000, slurm_partition='largemem'
    script:
        'scripts/annotate.py'

rule atac_model:
    input:
        cell_annotate = data_dir+'data/rna_cell_annot.csv',
        atac_anndata = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        umap_data = work_dir+'/data/atac_umap.csv',
        var_data = work_dir+'/data/atac_var_selected.csv',
        merged_atac_anndata = work_dir+'/atlas/03_filtered_anndata_atac.h5ad'
    params:
        samples=samples
    singularity:
        envs['atac']
    resources:
        runtime=2880, mem_mb=3000000, slurm_partition='largemem'
    threads:
        64
    script:
        work_dir+'/scripts/merge_atac.py' 

"""
# THIS IS EXPERIMENTAL AND NOT QUITE WORKING AND HARDCODED TO TEST EXAMPLE
rule atac_model:
    input:
        atac_anndata = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac.h5ad', 
            batch=batches,
            sample=samples
            )
    output:
        atac_umap = work_dir+'/data/atac_umap.csv',
        atac_var = work_dir+'/data/atac_var_selected.csv',
        
    singularity:
        envs['atac']
    threads:
        64
    resources:
        runtime=2880, disk_mb=500000, mem_mb=300000, gpu=4
    script:
        'scripts/atac_model.py'
"""

rule atac_annotate:
    input:
        atac_anndata = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac.h5ad', 
            zip,
            batch=batches,
            sample=samples
            ),
        umap_csv = work_dir+'/data/atac_umap.csv',
        var_csv = work_dir+'/data/atac_var_selected.csv',
        annot_csv = work_dir+'/data/rna_cell_annot.csv',
        metadata_table = metadata_table
    output:
        atac_anndata = work_dir+'/atlas/04_filtered_anndata_atac.h5ad',
        merged_atac_anndata = work_dir+'/atlas/05_annotated_anndata_atac.h5ad'
    params:
        samples=samples
    singularity:
        envs['atac']
    resources:
        runtime=2880, disk_mb=500000, mem_mb=300000
    script:
        'scripts/atac_annotate.py'

rule multiome_output:
    input:
        merged_atac_anndata = work_dir+'/atlas/05_annotated_anndata_atac.h5ad',
        merged_rna_anndata = work_dir+'/atlas/05_annotated_anndata_rna.h5ad'
    output:
        merged_multiome = work_dir+'/atlas/multiome_atlas.h5mu'
    singularity:
        envs['singlecell']
    script:
        'scripts/merge_muon.py'

rule export_celltypes:
    input:
        merged_multiome = work_dir+'/atlas/multiome_atlas.h5mu'
    output:
        celltype_atac = work_dir+'data/celltypes/{cell_type}/atac.h5ad',
        celltype_rna = work_dir+'data/celltypes/{cell_type}/rna.h5ad'
    params:
        cell_type = lambda wildcards, output: output[0].split('/')[-2]
    singularity:
        envs['singlecell']
    threads:
        8
    resources:
        runtime=120, mem_mb=300000
    script:
        'scripts/export_celltype.py'

rule DGE:
    input:
        rna_anndata = work_dir + '/atlas/05_annotated_anndata_rna.h5ad'
    output:
        output_DGE_data = work_dir + '/data/significant_genes/rna/rna_{cell_type}_{disease}_DGE.csv',
        output_figure = work_dir + 'figures/{cell_type}/rna_{cell_type}_{disease}_DAR.png',
        celltype_pseudobulk = work_dir+'/data/celltypes/{cell_type}/rna_{cell_type}_{disease}_pseudobulk.h5ad'
    params:
        disease_param = disease_param,
        control = control,
        disease = lambda wildcards, output: output[0].split("_")[-2],
        cell_type = lambda wildcards, output: output[0].split("_")[-3]
    singularity:
        envs['singlecell']
    threads:
        64
    resources:
        runtime=1440, disk_mb=200000, mem_mb=200000
    script:
        'scripts/rna_DGE.py'

rule DAR:
    input:
        atac_anndata = work_dir+'data/celltypes/{cell_type}/atac.h5ad'
    output:
        output_DAR_data = work_dir+'/data/significant_genes/atac/atac_{cell_type}_{disease}_DAR.csv',
        output_figure = work_dir+'/figures/{cell_type}/atac_{cell_type}_{disease}_DAR.png'
    params:
        disease_param = disease_param,
        control = control,
        disease = lambda wildcards, output: output[0].split("_")[-2],
        cell_type = lambda wildcards, output: output[0].split("_")[-3]
    singularity:
        envs['atac']
    threads:
        64
    resources:
        runtime=1440, disk_mb=200000, mem_mb=200000
    script:
        'scripts/atac_DAR.py'
