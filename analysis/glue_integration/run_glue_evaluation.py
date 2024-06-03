from evaluation import *
import anndata as ad
import pandas as pd
import scanpy as sc
import scglue
import networkx as nx
from itertools import chain

from configparser import ConfigParser, ExtendedInterpolation
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", help="Path to the config file.")
args = parser.parse_args()
config_path = Path(args.config)

config = ConfigParser(interpolation=ExtendedInterpolation())
config.read(config_path)

from matplotlib import rcParams
from matplotlib import pyplot as plt

scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)

data_dir = config['Paths']['data_dir']
model_name = config['Paths']['model_name']
annotation_file = config['Paths']['annotation_file']

HVGS = config['Features']['HVGS']
SELECTED_ATAC_FEATURES = config['Features']['SELECTED_ATAC_FEATURES']
model_path = "HVG_" + str(HVGS) + "_selected_atac_" + str(SELECTED_ATAC_FEATURES)

dataset = config['Paths']['dataset']

umap_file = data_dir + "/" + model_path + "/UMAP.pdf"
rna = ad.read_h5ad(data_dir + "/rna-pp.h5ad")
atac = ad.read_h5ad(data_dir + "/atac-pp.h5ad")

glue = scglue.models.load_model(data_dir + "/" + model_path + model_name)

rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)
combined = ad.concat([rna, atac], label="domain", keys=["RNA", "ATAC"])

annotation = pd.read_csv(annotation_file, delimiter='\t')

if dataset == "10XPBMC":
    annotation.index = annotation.rna_barcode.values
    if "CellinFile1#" in atac.obs_names[1]:
        annotation.index = "CellinFile1#" + annotation.rna_barcode.values
        if "CellinFile1#" not in rna.obs_names[1]:
            rna.obs_names = "CellinFile1#" + rna.obs_names
    atac.obs["rna_label"] = atac.obs.merge(annotation, left_index=True, right_index=True).final_label.values
    w = 0.8
if dataset == "Chen2019":
    annotation.index = annotation.barcode.values
    if "CellinFile1#" in atac.obs_names[1]:
        annotation.index = "CellinFile1#" + annotation.barcode.values
        if "CellinFile1#" not in rna.obs_names[1]:
            rna.obs_names = "CellinFile1#" + rna.obs_names
    atac.obs["rna_label"] = atac.obs.merge(annotation, left_index=True, right_index=True).final_label.values
    w = 0.4
    rna.obs["rna_label"]  = rna.obs["final_label"]



rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)
combined = ad.concat([rna, atac], label="domain", keys=["RNA", "ATAC"])
sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)

with plt.rc_context():
    sc.pl.umap(combined, color=["rna_label", "domain"], wspace=w, show=False)
    plt.savefig(umap_file, bbox_inches="tight")

rna = rna[find_indices(atac.obs_names, rna.obs_names),:]

# fracs = calc_domainAveraged_FOSCTTM(rna.obsm["X_glue"], atac.obsm["X_glue"])
# df_fracs = pd.DataFrame({
#         'barcode': rna.obs_names,
#         'FOSCTTM': fracs
#     })
# df_fracs.to_csv(data_dir + "/" + model_path + '/FOSCTTM.csv', 
#                  index=False)

fracs, a, b = calc_domainAveraged_FOSCTTM2(rna.obsm["X_glue"], atac.obsm["X_glue"], dist_name="cosine", labels=rna.obs['rna_label'])
df_fracs = pd.DataFrame({
        'barcode': rna.obs_names,
        'FOSCTTM': fracs,
        'FOSCTTM_same':a,
        'FOSCTTM_diff':b
    })
df_fracs.to_csv(data_dir + "/" + model_path + '/FOSCTTM_cosine.csv', 
                 index=False)