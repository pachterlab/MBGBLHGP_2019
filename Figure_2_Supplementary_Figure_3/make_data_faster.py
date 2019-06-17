#!/usr/bin/env python
# coding: utf-8

import argparse
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import time
import copy
from sklearn.preprocessing import LabelEncoder
from scipy import sparse
import scipy
import anndata
from matplotlib.pyplot import figure
from sklearn.decomposition import TruncatedSVD
import sklearn
import anndata
import time
from openTSNE import TSNE
import openTSNE
from openTSNE.callbacks import ErrorLogger
import datetime
from anndata import AnnData
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.pairwise import manhattan_distances
import matplotlib.pyplot as plt
import pickle

def make_data_faster(dataset_shortname):
	k_folder = '/home/single_cell_analysis/kallisto_out_single_bustools_dev/kallisto_' + dataset_shortname
	if dataset_shortname in ["pbmc_1k_v3", "pbmc_10k_v3", "neuron_10k_v3"]:
	     dataset_shortname = dataset_shortname.split("_")[0] + dataset_shortname.split("_")[1] + "_" + dataset_shortname.split("_")[2]
	c_folder = '/home/single_cell_analysis/cellranger_out/cellranger3_' + dataset_shortname +'_out/outs/filtered_feature_bc_matrix'
	c_raw_folder = '/home/single_cell_analysis/cellranger_out/cellranger3_' + dataset_shortname +'_out/outs/raw_feature_bc_matrix'

	c_raw = anndata.AnnData(scipy.io.mmread(os.path.join(c_raw_folder,'matrix.mtx.gz')).tocsr().T)
	c_barcodes = pd.read_csv(os.path.join(c_raw_folder,'barcodes.tsv.gz'), index_col = 0, header = None, names = ['barcode'])
	c_barcodes.index = c_barcodes.index.str.slice(0,16,1)
	c_raw.obs = c_barcodes
	c_raw.var = pd.read_csv(os.path.join(c_raw_folder,'features.tsv.gz'), header = None, index_col = 0, names =['ensembl_id', 'gene_name', 'kind'], sep = '\t')
	print('Loaded c raw mtx:',c_raw.X.shape)

	del c_barcodes

	# load c filtered matrix
	c = anndata.AnnData(scipy.io.mmread(os.path.join(c_folder,'matrix.mtx.gz')).tocsr().T)
	c_barcodes = pd.read_csv(os.path.join(c_folder,'barcodes.tsv.gz'), index_col = 0, header = None, names = ['barcode'])
	c_barcodes.index = c_barcodes.index.str.slice(0,16,1)
	c.obs = c_barcodes
	c.var = pd.read_csv(os.path.join(c_folder,'features.tsv.gz'), header = None, index_col = 0, names =['ensembl_id', 'gene_name', 'kind'], sep = '\t')
	print('Loaded c filtered mtx:',c.X.shape)

	del c_barcodes


	## load kallisto raw matrix
	k_raw = anndata.AnnData(scipy.io.mmread(os.path.join(k_folder,'genes.mtx')).tocsr())
	k_raw.obs= pd.read_csv(os.path.join(k_folder,'genes.barcodes.txt'), index_col = 0, header = None, names = ['barcode'])
	k_raw.var = pd.read_csv(os.path.join(k_folder,'genes.genes.txt'), header = None, index_col = 0, names =['ensembl_id'], sep = '\t')
	print('Loaded k raw mtx:',k_raw.X.shape)


	# truncdates the ensembl version number off the kallisto labels
	k_raw.var['full_emsembl_id'] = k_raw.var.index
	k_raw.var.index = k_raw.var['full_emsembl_id'].str.slice(0,18)


	if dataset_shortname in ['hgmm1k_v2', 'hgmm1k_v3', 'hgmm10k_v3']:
	    k_raw.var.index = k_raw.var['full_emsembl_id']

	    # do this as late as possible
	k = k_raw[c.obs.index.values]
	print('Loaded k filtered mtx:', k.X.shape)

	c_raw.obs['counts'] = c_raw.X.sum(1)
	c_raw.obs['ngenes'] = np.array((c_raw.X > 0).sum(1))
	c_raw = c_raw[c_raw.obs['counts'] > 0]
	c_raw.layers['log1p'] = np.log1p(c_raw.X)
	c_raw.obs['log10counts']= np.log10(c_raw.obs['counts'])
	print('Cell Ranger raw:', c_raw.shape)


	# count UMIs, genes, log transform raw kallisto barcodes   
	# first remove kallisto barcodes with 0 gene counts

	k_raw.obs['counts'] = k_raw.X.sum(1)
	k_raw.obs['ngenes'] = np.array((k_raw.X > 0).sum(1))
	k_raw = k_raw[k_raw.obs['counts'] > 0]
	k_raw.layers['log1p'] = np.log1p(k_raw.X)
	k_raw.obs['log10counts'] = np.log10(k_raw.obs['counts'])
	print('kallisto raw:', k_raw.shape)

	c.obs['counts'] = c.X.sum(1)
	c.obs['ngenes'] = np.array((c.X > 0).sum(1))
	c = c[c.obs['counts'] > 0]
	c.layers['log1p'] = np.log1p(c.X)
	c.obs['log10counts']= np.log10(c.obs['counts'])
	print('Cell Ranger filtered:', c.shape)


	# count UMIs, genes, log transform filtered kallisto barcodes   
	# first remove kallisto barcodes with 0 gene counts

	k.obs['counts'] = k.X.sum(1)
	k.obs['ngenes'] = np.array((k.X > 0).sum(1))
	k = k[k.obs['counts'] > 0]
	k.layers['log1p'] = np.log1p(k.X)
	k.obs['log10counts'] = np.log10(k.obs['counts'])
	print('kallisto filtered:', k.shape)

	joint_obs = k_raw.obs.join(c_raw.obs,how = 'outer', lsuffix='-kallisto', rsuffix='-tenx')
	joint_obs = joint_obs.fillna(0)
	print('Total barcodes seen')
	print(len(joint_obs))

	# barcodes seen by both
	common_obs = k_raw.obs.join(c_raw.obs,how = 'inner', lsuffix='-kallisto', rsuffix='-tenx')
	print('Barcodes seen by both')
	print(len(common_obs))

	kobs = k_raw.obs.join(c_raw.obs,how = 'left', lsuffix='-kallisto', rsuffix='-tenx')
	kobs = kobs.sort_values(by=['counts-kallisto'], ascending = False)
	print('Barcodes seen by kallisto missed by Cell Ranger')
	print(len(joint_obs) - len(kobs))


	# just Cell Ranger observations
	tobs = c_raw.obs.copy()
	tobs = tobs.sort_values('counts', ascending = False)
	print('Barcodes seen by Cell Ranger missed by kallisto')
	print(len(joint_obs) - len(tobs))

	# ## Compute correlations between kallisto and Cell Ranger
	# handy and fast function for computing correlation on sparse matrices
	def sparse_M_std(X):
	    n = X.shape[1]
	    return np.sqrt(n * X.multiply(X).sum(1) - np.multiply(X.sum(1), X.sum(1)))

	def sparse_M_corr(X,Y):
	    X_std = sparse_M_std(X)
	    Y_std = sparse_M_std(Y)
	    XY_std = np.multiply(X_std, Y_std)

	    n = X.shape[1]
	    XY_cov = n* X.multiply(Y).sum(1) - np.multiply(X.sum(1), Y.sum(1))
	    R = np.divide(XY_cov, XY_std)
	    return np.squeeze(np.asarray(R))

	raw_counts_correlation = sparse_M_corr(k_raw[common_obs.index].layers['log1p'],c_raw[common_obs.index].layers['log1p'])
	filtered_counts_correlation = sparse_M_corr(k_raw[c.obs.index].layers['log1p'],c_raw[c.obs.index].layers['log1p'])
	print('Correlations computed!')

	tsvd = TruncatedSVD(n_components=10)
	TSVD = tsvd.fit_transform(k.layers['log1p'])
	k.obsm['TSVD'] = TSVD
	k.obsm['TSVD']
	print('TSVD variance ratios:\n', list(tsvd.explained_variance_ratio_))
	print(datetime.datetime.now())


	tsvd = TruncatedSVD(n_components=10)
	TSVD = tsvd.fit_transform(c.layers['log1p'])
	c.obsm['TSVD'] = TSVD
	c.obsm['TSVD']
	print('TSVD variance ratios:\n', list(tsvd.explained_variance_ratio_))
	print(datetime.datetime.now())


	print('Calculating L1 distances...')

	# taking manhattan distance between matrices
	dnck = manhattan_distances(c.layers['log1p'], k.layers['log1p'])
	dnkk = manhattan_distances(k.layers['log1p'], k.layers['log1p'])
	print(datetime.datetime.now())

	# nkc are the kallisto-cellranger distances 
	nck = np.diagonal(dnck)

	# ncc are the kallisto-kallisto distances
	nkk = []
	for row in dnkk:
	    val = np.partition(row, 1)[1]
	    nkk.append(val)
	print('L1 distances done!')
	print(datetime.datetime.now())


	print('Doing t-SNE')
	print(datetime.datetime.now())
	tsne = TSNE(perplexity=30, metric="euclidean", callbacks=openTSNE.callbacks.ErrorLogger(),n_jobs=8, random_state=42, n_iter=750 )
	k.obsm['TSNE10'] = tsne.fit(k.obsm['TSVD'])
	print('kallisto TSNE-10 done.')
	print(datetime.datetime.now())


	# Perform TSNE on top 10 truncated SVD components of Cell Ranger filtered matrix

	print('Doing t-SNE on top 10 PC for Cell Ranger')
	# 
	print(datetime.datetime.now())
	tsne = TSNE(perplexity=30, metric="euclidean", callbacks=openTSNE.callbacks.ErrorLogger(),n_jobs=8, random_state=42, n_iter=750 )
	c.obsm['TSNE10'] = tsne.fit(c.obsm['TSVD'])
	print('Cell Ranger TSNE-10 done.')
	print(datetime.datetime.now())


	c_raw.write(os.path.join("./write_data/" + dataset_shortname + '_tenx_raw.h5ad'))
	k_raw.write(os.path.join("./write_data/" + dataset_shortname + '_kallisto_raw.h5ad'))
	k.write(os.path.join("./write_data/" + dataset_shortname + '_kallisto.h5ad'))
	c.write(os.path.join("./write_data/" + dataset_shortname + '_tenx.h5ad'))


	with open(os.path.join("./write_data/" + dataset_shortname + '_kobs.pkl'), 'wb') as handle:
	    pickle.dump(kobs, handle, protocol=pickle.HIGHEST_PROTOCOL)
	    
	with open(os.path.join("./write_data/" + dataset_shortname + '_tobs.pkl'), 'wb') as handle:
	    pickle.dump(tobs, handle, protocol=pickle.HIGHEST_PROTOCOL)
	    
	with open(os.path.join("./write_data/" + dataset_shortname + '_common_obs.pkl'), 'wb') as handle:
	    pickle.dump(common_obs, handle, protocol=pickle.HIGHEST_PROTOCOL)

	with open(os.path.join("./write_data/" + dataset_shortname + '_joint_obs.pkl'), 'wb') as handle:
	    pickle.dump(joint_obs, handle, protocol=pickle.HIGHEST_PROTOCOL)
	    
	with open(os.path.join("./write_data/" + dataset_shortname + '_nkk.pkl'), 'wb') as handle:
	    pickle.dump(nkk, handle, protocol=pickle.HIGHEST_PROTOCOL)
	    
	with open(os.path.join("./write_data/" + dataset_shortname + '_nck.pkl'), 'wb') as handle:
	    pickle.dump(nck, handle, protocol=pickle.HIGHEST_PROTOCOL)
	    
	with open(os.path.join("./write_data/" + dataset_shortname + '_raw_counts_correlation.pkl'), 'wb') as handle:
	    pickle.dump(raw_counts_correlation, handle, protocol=pickle.HIGHEST_PROTOCOL)
	    
	with open(os.path.join("./write_data/" + dataset_shortname + '_filtered_counts_correlation.pkl'), 'wb') as handle:
	    pickle.dump(filtered_counts_correlation, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="processes all data")
    parser.add_argument("--ds", help="dataset name (i.e. SRR8599150_v2)")
    

    args = parser.parse_args()
    print("Loading files..")
    make_data_faster(str(args.ds))

    print("Done")
