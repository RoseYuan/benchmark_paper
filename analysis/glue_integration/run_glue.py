from itertools import chain

import anndata as ad
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns

from configparser import ConfigParser, ExtendedInterpolation
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", help="Path to the config file.")
args = parser.parse_args()
config_path = Path(args.config)

config = ConfigParser(interpolation=ExtendedInterpolation())
config.read(config_path)

HVGS = config['Features']['HVGS']
SELECTED_ATAC_FEATURES = config['Features']['SELECTED_ATAC_FEATURES']

# dataset = config['Paths']['dataset']
# method = config['Paths']['method']
# out_path = config['Paths']['out_path']
model_path = "HVG_" + str(HVGS) + "_selected_atac_" + str(SELECTED_ATAC_FEATURES)
model_name = config['Paths']['model_name']
data_dir = config['Paths']['data_dir']

rna = ad.read_h5ad(data_dir + "/rna-pp.h5ad")
atac = ad.read_h5ad(data_dir + "/atac-pp.h5ad")

if SELECTED_ATAC_FEATURES:
    atac = atac[:,atac.uns['var_features'].index]
    guidance = nx.read_graphml(data_dir + "/guidance_var.graphml.gz")
else:
    guidance = nx.read_graphml(data_dir + "/guidance.graphml.gz")


scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=HVGS,
    use_layer="counts", use_rep="X_pca"
)

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=HVGS,
    use_rep="embedding"
)

if HVGS:
    guidance = guidance.subgraph(chain(rna.var.query("highly_variable").index, atac.var.query("highly_variable").index)).copy()


glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, guidance,
    fit_kws={"directory": data_dir + "/" + model_path}
)


glue.save(data_dir + "/" + model_path + model_name)


# glue = scglue.models.load_model("/home/siluo/public/SiyuanLuo/projects/rebuttal/integration/outputs/10XPBMC/ArchR_peaks/all_selected_features_4.22/glue.4.22.dill")

# dx = scglue.models.integration_consistency(
#     glue, {"rna": rna, "atac": atac}, guidance, n_meta=100
# )

# print()
