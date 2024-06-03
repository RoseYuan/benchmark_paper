import anndata as ad
import networkx as nx
import scanpy as sc
import scglue

from configparser import ConfigParser, ExtendedInterpolation
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", help="Path to the config file.")
args = parser.parse_args()
config_path = Path(args.config)

config = ConfigParser(interpolation=ExtendedInterpolation())
config.read(config_path)

###############################################
print("Load data.....")

# dataset = config['Paths']['dataset']
method = config['Paths']['method']
# out_path = config['Paths']['out_path']
# rna_path = out_path + dataset + "/rna.h5ad"
# atac_path = out_path + dataset + "/" + method + "/loaded_obj.h5ad"
gtf_path = config['Paths']['gtf_path']
gtf_by = config['Paths']['gtf_by']
data_dir = config['Paths']['data_dir']
atac_path = data_dir + "/loaded_obj.h5ad"
rna_path = data_dir + "/../rna.h5ad"

rna = ad.read_h5ad(rna_path)
atac = ad.read_h5ad(atac_path)

################################################
print("Annotate features.....")

# annotate genes
scglue.data.get_gene_annotation(
    rna, gtf=gtf_path,
    gtf_by=gtf_by
)

rna = rna[:, rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].dropna().index]
rna.var['chromStart'] = rna.var['chromStart'].astype(int)
rna.var['chromEnd'] = rna.var['chromEnd'].astype(int)
rna.var['chrom'] = "chr" + rna.var['chrom'].astype(str)

split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)

#################################################
print("Get variable features.....")

if "ArchR" in method:
    atac.uns['var_features'].index = atac.uns['var_features']['seqnames'] + \
        ':' + atac.uns['var_features']['start'].astype('Int64').astype(str) + '-' \
        + atac.uns['var_features']['end'].astype('Int64').astype(str)

# subsetting
atac_var_features = atac[:,atac.uns['var_features'].index]

#################################################
# print("Construct graph using all genes and all atac features.....")
# guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
# scglue.graph.check_graph(guidance, [rna, atac])

print("Construct graph using all genes and selected atac features.....")
rna.var.highly_variable = rna.var.highly_variable.astype('bool')
guidance_var = scglue.genomics.rna_anchored_guidance_graph(rna, atac_var_features)
scglue.graph.check_graph(guidance_var, [rna, atac_var_features])

rna.write(data_dir + "/rna-pp.h5ad", compression="gzip")
atac_var_features.write(data_dir + "/atac-pp.h5ad", compression="gzip")
# nx.write_graphml(guidance, data_dir + "/guidance.graphml.gz")
nx.write_graphml(guidance_var, data_dir + "/guidance_var.graphml.gz")