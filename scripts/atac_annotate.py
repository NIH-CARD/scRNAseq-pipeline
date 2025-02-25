import snapatac2 as snap
import pandas as pd

# Import the samples
samples = snakemake.params.samples
atac_anndata = snakemake.input.atac_anndata

# Read in snapATAC2 datasets into a list of anndata objects in read only
list_of_anndata = [(samples[i], snap.read(atac_anndata[i])) for i in range(len(atac_anndata))]

# Create the AnnDataSet from the snapATAC2 datasets
anndataset = snap.AnnDataSet(
    adatas=list_of_anndata,
    filename=atac_anndata
)

#
metadata_dict = {
    'batch' : 'Use_batch',
    'sex': 'Sex',
    'age': 'Age',
    'brain_bank': 'Brain_bank',
    'short diagnosis': 'Short Diagnosis',
    'Primary Diagnosis': 'Primary Diagnosis'
}

metadata_df = pd.read_csv(snakemake.input.metadata_table)

for new_obs, metadata_obs in metadata_dict.items():
    print(new_obs)
    anndataset.obs[new_obs] = [metadata_df.loc[metadata_df['Sample_ID'] == x][metadata_obs].iloc[0] for x in anndataset.obs['sample']]

# Read in the UMAP previously calculated
umap_df = pd.read_csv(snakemake.input.umap_csv)
anndataset.obsm['X_umap'] = umap_df[['umap_x', 'umap_y']].to_numpy()

# Input the total count and selected 
var_df = pd.read_csv(snakemake.input.var_csv)
anndataset.var['count'] = var_df['count'].to_numpy()
anndataset.var['selected'] = var_df['selected'].to_numpy()

# Import the RNA annotation
rna_annot = pd.read_csv(snakemake.input.annot_csv)
anndataset.obs['cell_type'] = rna_annot['cell_type'].to_list()

atac = anndataset.to_adata()

atac.write_h5ad(snakemake.output.merged_atac_anndata, compression='gzip')