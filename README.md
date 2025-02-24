# multiome-pipeline

## Single-Nuclei Multiome Atlas Pipeline

This pipeline is built off of the Substantia Nigra single-nuclei atlas pipeline here (https://github.com/adamcatching/SN_atlas) and is being updated to be more interchangable between projects. Starting off the scripts and Snakemake pipeline work for the Substantia Nigra dataset.

The multiome pipeline was built to take several batches of human brain single-nuclei sequencing samples and process them into a multiome atlas object for further analysis. There are a currently number of hardcoded parameters such as the raw data paths and which cell type markers to be used for annotation. These should be modified for running this pipeline on a separate atlas, future plans are for making a more universal multiome pipeline with a single config file (see https://github.com/adamcatching/multiome-pipeline). The modules of Scanpy (https://github.com/scverse/scanpy), SCVI (https://github.com/scverse/scvi-tools), and snapATAC2 (https://github.com/kaizhang/SnapATAC2) are utilized heavily to produce a multiome atlas with minimal batch effects.

Required inputs:
- Metadata file in .csv format
- Cell-typing table with marker genes, in .csv format
- Directory of the CellRanger output

Outputs:
- RNA and ATAC Multiome atlas object (multiome_atlas.h5ad)
- List of differentially expressed genes and accessible regions (data/significant_genes/(rna or atac)/(celltype)_(disease)_(DGE or DAR).csv
- Preprocessed and QC-filtered AnnData objects for each sample

Current version:
- Uses Singularity images for reproducible runs
- Snakemake runs all steps until all output files are created
- Genes used for cell-typing are hard-coded into the Annotate rule
- Both RNA and ATAC processing has to be done to run this pipeline 
- Requires a metadata .csv file and all values within are saved in the atlas.
- Data needs to be stored in a specific heirarchy
- Cellbender needs to be run after CellRanger; _once ambient RNA is corrected in future CellRanger-ARC this won't be necessary_.
- Differential Gene Expression and Differential Accessibility of Regions analysis are done in separate notebooks; _these may be integrated into a rule with differential parameters specified in the config file_.

## Pipeline

![screenshot](images/multiome_pipeline.png)

Once set up, this complete pipeline can be run by simply typing '''bash snakemake.sh ''' in terminal in an HPC running Slurm. This is a work in progress and has not been tested on other devices. 

### RNA processing

Transcriptomic data from CellRanger-ARC-2.0 ('''cellbender_gex_counts_filtered.h5''') is read in and processed with Scanpy. QC metrics of percent mitochondria/ribosomal RNA, doublet probability, and cell cycle.

### RNA QC

Parameters from the processing step are used to filter the cells from each samples based on percent mitochondrial transcripts, probability of being a doublet, and the minimum number of genes observed per cell.

### Individual RNA sample merging into atlas

Each individual RNA AnnData object are merged into a single QC-filtered object for downstream analysis.

### ATAC processing

ATAC fragment data is converted into an AnnData object with bins used as the measured variable in each cell. One object is created for each sample.

### ATAC QC

Cells in each sample's ATAC object are filtered for a minimum number of bins per cell. 

### Filtering RNA and ATAC data 

Each sample's QC-filtered RNA and ATAC AnnData objects are filtered for the same cells observed in both samples. Final AnnData objects are saved with a '''03_''' prefix.

### RNA modeling

Filtered RNA samples are merged into an atlas and multidimensional scaling is performed. A copy of the atlas is made with mitochondiral and ribosomal transcripts removed and only the most variable genes kept. SCVI is used to model the embed dimensions of the atlas, with batch correction, followed by KNN, leiden clustering, and UMAP scaling.

### Cell-typing

Cell types of the modeled and clustered RNA atlas are estimated using over-representation analysis and a currated list of cell gene markers.

### ATAC modeling

Using snapATAC2 for only read/write, both create a AnnDataSet object to run batch-corrected spectral analysis and scaling; resulting in a leiden-clustered UMAP AnnData object

### Merging to one multiome object
Both atlases are merged into a single muon AnnData object for portability.
