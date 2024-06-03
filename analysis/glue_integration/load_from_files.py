import os
import pandas as pd
import numpy as np

import scanpy as sc

from configparser import ConfigParser, ExtendedInterpolation
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", help="Path to the config file.")
args = parser.parse_args()

config_path = Path(args.config)
config = ConfigParser(interpolation=ExtendedInterpolation())
config.read(config_path)


# dataset = config['Paths']['dataset']
method = config['Paths']['method']
# out_path = config['Paths']['out_path']
# input_path = config['Paths']['input_path']
data_dir = config['Paths']['data_dir']
embedding_file = config['Paths']['embedding_file']
# data_dir =  str(out_path)  + dataset + "/" + method

from scipy.io import mmread
counts = mmread(data_dir + "/counts.mtx")

cells = pd.read_csv(data_dir + '/cells.csv', index_col=0).iloc[:, 0]
features = pd.read_csv(data_dir + '/features.csv', index_col=0)
var_features = pd.read_csv(data_dir + '/LSIFeatures.csv', index_col=0)

if method == "ArchR_tiles":
    features['end'] = features['start'] + 500
    var_features['end'] = var_features['start'] + 500
if 'ArchR' in method:
    features.index = features['seqnames'] + ':' + features['start'].astype('Int64').astype(str) + '-' + features['end'].astype('Int64').astype(str)
    var_features.index = var_features['seqnames'] + ':' + var_features['start'].astype('Int64').astype(str) + '-' + var_features['end'].astype('Int64').astype(str)
if ('Signac' in method) or (method == "aggregation"):
    var_features.index = var_features.x.values

# counts and indices
ad = sc.AnnData(counts.T)
ad.obs_names = cells
ad.var_names = features.index

# feature meta data
for col in features.columns:
    ad.var[col] = features[col]
ad.X = ad.X.tocsr()

# low dimensional embeddings
ad.obsm['embedding'] = pd.read_csv(embedding_file, index_col=0, delimiter='\t', header=None).loc[ad.obs_names, : ].values

# cell meta data
cell_meta = pd.read_csv(data_dir + '/cell_metadata.csv', index_col=0).loc[ad.obs_names, : ]
for col in cell_meta.columns:
    ad.obs[col] = cell_meta[col].values

# variable features
ad.uns['var_features'] = var_features

out_file = data_dir + '/loaded_obj.h5ad'
ad.write_h5ad(out_file)