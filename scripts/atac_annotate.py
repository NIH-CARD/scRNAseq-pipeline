import scanpy as sc
import pandas as pd

adata = sc.read_h5ad(snakemake.input.merged_atac_anndata)

# Import the RNA annotation
rna_annot = pd.read_csv(snakemake.input.annot_csv)
# Convert to a dictionary
celltype_bc = dict(zip(rna_annot['atlas_identifier', 'cell_type']))

# Assign cell types from RNA data
adata.obs['cell_type'] = [celltype_bc[x] for x in adata.obs['atlas_identifier']]

adata.write_h5ad(snakemake.output.merged_atac_anndata, compression='gzip')