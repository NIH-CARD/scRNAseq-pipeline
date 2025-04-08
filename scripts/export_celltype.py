import muon as mu
import scanpy as sc

# Read in data
mdata = mu.read_h5mu(snakemake.input.merged_multiome)

# Split out only the celltype AnnData, and reduce the size of the object
mdata = mdata[mdata.obs['rna:cell_type'] == snakefile.params.cell_type].copy()

# Split out the RNA AnnData
rna = mdata['rna'].copy()

# Split out the ATAC AnnData
atac = mdata['atac'].copy()

# Write the rna data out
sc.write_h5ad(snakefile.output.celltype_rna, compression='gzip')

# Write the atac data out
sc.write_h5ad(snakefile.output.cell_type_atac, compression='gzip')
