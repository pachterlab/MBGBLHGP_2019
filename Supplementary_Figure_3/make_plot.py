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

def make_plot(dataset_shortname):
    # alldata = anndata.read_h5ad(os.path.join('./write_data/' + dataset_shortname + '_alldata.h5ad'))
    tenx_raw = anndata.read_h5ad(os.path.join('./write_data/' + dataset_shortname + '_tenx_raw.h5ad'))
    kallisto_raw = anndata.read_h5ad(os.path.join('./write_data/' + dataset_shortname + '_kallisto_raw.h5ad'))
    kallisto = anndata.read_h5ad(os.path.join('./write_data/' + dataset_shortname + '_kallisto.h5ad'))
    tenx = anndata.read_h5ad(os.path.join('./write_data/' + dataset_shortname + '_tenx.h5ad'))
    
    # kallisto_copy = anndata.read_h5ad(os.path.join('./write_data/' + dataset_shortname + '_kallisto_copy.h5ad'))
    # tenx_copy = anndata.read_h5ad(os.path.join('./write_data/' + dataset_shortname + '_tenx_copy.h5ad'))
    
    with open(os.path.join("./write_data/" + dataset_shortname + '_kobs.pkl'), 'rb') as handle:
        kobs = pickle.load(handle)
        
    with open(os.path.join("./write_data/" + dataset_shortname + '_tobs.pkl'), 'rb') as handle:
        tobs = pickle.load(handle)

    with open(os.path.join("./write_data/" + dataset_shortname + '_common_obs.pkl'), 'rb') as handle:
        common_obs = pickle.load(handle)

    with open(os.path.join("./write_data/" + dataset_shortname + '_joint_obs.pkl'), 'rb') as handle:
        joint_obs = pickle.load(handle)

    with open(os.path.join("./write_data/" + dataset_shortname + '_nkk.pkl'), 'rb') as handle:
        nkk = pickle.load(handle)
        
    with open(os.path.join("./write_data/" + dataset_shortname + '_nck.pkl'), 'rb') as handle:
        nck = pickle.load(handle)
        
    with open(os.path.join("./write_data/" + dataset_shortname + '_raw_counts_correlation.pkl'), 'rb') as handle:
        raw_counts_correlation = pickle.load(handle)
        
    with open(os.path.join("./write_data/" + dataset_shortname + '_filtered_counts_correlation.pkl'), 'rb') as handle:
        filtered_counts_correlation = pickle.load(handle)
    
    
    def lighten_color(color, amount=0.5):
        """
        Lightens the given color by multiplying (1-luminosity) by the given amount.
        Input can be matplotlib color string, hex string, or RGB tuple.

        Examples:
        >> lighten_color('g', 0.3)
        >> lighten_color('#F034A3', 0.6)
        >> lighten_color((.3,.55,.1), 0.5)
        """
        import matplotlib.colors as mc
        import colorsys
        try:
            c = mc.cnames[color]
        except:
            c = color
        c = colorsys.rgb_to_hls(*mc.to_rgb(c))
        return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


    # In[24]:


    # define some values we use to keep the plots consistent
    fsize= 15
    kallisto_color = '#FF7F0E'
    tenx_color = '#1F77B4'
    dotsize = 10
    xmax = 1e5
    gridalpha = 0.2

    # Make the figure with a 4x4 grid
    # This is necessary so panel G can have left and right sides
    fig = plt.figure(figsize=(21,25))
    ax_a = plt.subplot2grid((4, 4), (0, 0), colspan=2)
    ax_b = plt.subplot2grid((4, 4), (1, 0), colspan=2)
    ax_c = plt.subplot2grid((4, 4), (2, 0), colspan=2)
    ax_d = plt.subplot2grid((4, 4), (3, 0), colspan=2)
    ax_e = plt.subplot2grid((4, 4), (0, 2), colspan=2)
    ax_f = plt.subplot2grid((4, 4), (1, 2), colspan=2)
    ax_g = plt.subplot2grid((4, 4), (2, 2), colspan=2)
    #ax_g_left = plt.subplot2grid((4, 4), (2, 2), colspan=1)
    #ax_g_right = plt.subplot2grid((4, 4), (2, 3), colspan=1)
    ax_h = plt.subplot2grid((4, 4), (3, 2), colspan=2)

    # Now we make the plots

    # A ###  KNEE PLOT ##########
    loc = ax_a

    tenx_ranked_umi = np.sort(np.array(tenx_raw.X.sum(1)), axis=None)[::-1]
    loc.plot( tenx_ranked_umi, np.arange(len(tenx_ranked_umi)),
             c = lighten_color(tenx_color, 0.5), linewidth=2, alpha=1)
    kallisto_ranked_umi = np.sort(np.array(kallisto_raw.X.sum(1)), axis=None)[::-1]
    loc.plot( kallisto_ranked_umi, np.arange(len(kallisto_ranked_umi)),  
             c = lighten_color(kallisto_color, 0.5), linewidth=2, alpha=1)
    tenx_ranked_umi = np.sort(np.array(tenx_raw.X.sum(1)), axis=None)[::-1]
    loc.plot( tenx_ranked_umi[0:np.shape(kallisto.X)[0]], np.arange(np.shape(kallisto.X)[0]),
             color =tenx_color, linewidth=2, label = 'Cell Ranger ', alpha=1)
    kallisto_ranked_umi = np.sort(np.array(kallisto_raw.X.sum(1)), axis=None)[::-1]
    loc.plot( kallisto_ranked_umi[0:np.shape(kallisto.X)[0]], np.arange(np.shape(kallisto.X)[0]),  color =kallisto_color, linewidth=2, label = 'kallisto', alpha=1)

    loc.set_xscale('log')
    loc.set_xlim(1,xmax) 
    loc.set_yscale("log", nonposy='clip')
    loc.set_xlabel('kallisto UMI counts',fontsize=fsize)
    loc.set_ylabel('Cumulative number of barcodes',fontsize=fsize)
    loc.set_title('',loc='center')
    loc.set_title('A', fontweight='bold', fontsize = fsize, loc = 'left' )
    loc.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = gridalpha)
    loc.axhline(y=np.shape(kallisto.X)[0],linewidth=2, color='black', linestyle='--')
    loc.axvline(x=kallisto_ranked_umi[np.shape(kallisto.X)[0]],linewidth=2, color='black', linestyle='--')

    handles, labels = loc.get_legend_handles_labels()
    loc.legend(handles[::-1], labels[::-1])

    # B ###  BARCODE RATIOS PLOT ##########
    loc = ax_b

    loc.plot(np.geomspace(1,10e5,100),np.geomspace(1,10e5,100),'gray',linewidth=1 ) # identity line
    loc.scatter(joint_obs['counts-kallisto'].values,joint_obs['counts-tenx'].values,   
                color ='lightgray', s=dotsize, alpha=0.3, edgecolors = 'none', label = 'Discarded barcodes')
    loc.scatter(kallisto.obs['counts'].values,tenx.obs['counts'].values,   
                color ='black', s=dotsize, alpha=0.3, edgecolors = 'none', label = 'Retained barcodes')

    # loc.set_aspect('equal')
    loc.set_xscale('log')
    loc.set_yscale("log", nonposy='clip')
    loc.set_xlabel('kallisto UMI counts',fontsize=fsize)
    loc.set_ylabel('Cell Ranger UMI counts',fontsize=fsize)
    loc.set_title('',loc='center')
    loc.set_xlim(1,xmax) 
    loc.set_ylim(1,xmax) 
    loc.set_title('B', fontweight='bold', fontsize = fsize, loc = 'left' )
    loc.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = gridalpha)
    handles, labels = loc.get_legend_handles_labels()
    loc.legend(handles[::-1], labels[::-1])

    # C ### GENE VS UMI PLOT ####
    loc = ax_c

    # discarded barcodes in lighter colors
    loc.scatter(kobs['counts-tenx'],kobs['ngenes-tenx'], s =dotsize, alpha=0.5, 
                c = lighten_color(tenx_color, 0.3), edgecolors = 'none', label  = 'Cell Ranger discarded barcodes' )
    loc.scatter(kobs['counts-kallisto'],kobs['ngenes-kallisto'], s = dotsize, alpha=0.5, 
                c = lighten_color(kallisto_color, 0.3), edgecolors = 'none', label  = 'kallisto discarded barcodes' )

    # retained barcodes in darker color
    loc.scatter(tenx.obs['counts'],tenx.obs['ngenes'], s =dotsize, alpha=0.5,
                c = tenx_color, edgecolors = 'none', label  = 'Cell Ranger retained barcodes' ) 
    loc.scatter(kallisto.obs['counts'],kallisto.obs['ngenes'], s = dotsize, alpha=0.5,
                c = kallisto_color, edgecolors = 'none', label  = 'kallisto retained barcodes' )

    loc.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = 0.5)
    loc.set_xscale('log')
    loc.set_xlim(1,xmax) 
    loc.set_yscale("log", nonposy='clip')
    loc.set_ylabel('Genes detected',fontsize=fsize)
    loc.set_xlabel('kallisto UMI counts',fontsize=fsize)
    loc.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = gridalpha)
    loc.set_title('C', fontweight='bold', fontsize = fsize, loc = 'left' )
    handles, labels = loc.get_legend_handles_labels()
    loc.legend(handles[::-1], labels[::-1])

    # D ### CORRELATION PLOT ####
    loc = ax_d

    # discarded barcodes in gray
    loc.scatter(x = common_obs['counts-kallisto'], y = raw_counts_correlation, s =dotsize, c = 'lightgray',
                alpha=0.3, edgecolors = 'none' , label = 'Discarded barcodes')
    # retained barcodes in black
    loc.scatter(x = kallisto.obs['counts'], y = filtered_counts_correlation, s =dotsize, c = 'black', 
                alpha=0.3, edgecolors = 'none', label = 'Retained barcodes' )

    loc.set_xscale('log')
    loc.set_xlim(1,xmax) 
    loc.set_ylim(0,1) 
    loc.set_title('D', fontweight='bold', fontsize = fsize, loc = 'left')
    loc.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = gridalpha)
    loc.set_xlabel('kallisto UMI counts', fontsize = fsize)
    loc.set_ylabel('Pearson Correlation', fontsize = fsize)
    handles, labels = loc.get_legend_handles_labels()
    loc.legend(handles[::-1], labels[::-1])

    # E ### BARCODE DISTANCE PLOT ####
    loc = ax_e

    # to choose the bins we pick the ones that cover the largest interval as np.histogram produces them
    # for nck, ncc, and the concatenation of ncc, nck
    hist, concat_bins = np.histogram(np.concatenate((nkk,nck)), bins='auto')
    hist, ck_bins =  np.histogram(nck, bins='auto')
    hist, cc_bins =  np.histogram(nkk, bins='auto')
    best_bins = max([ck_bins,concat_bins,cc_bins], key=len)

    loc.hist(x=nck, bins=best_bins, alpha=0.5, color = tenx_color, label="Closest Cell Ranger barcode")
    loc.hist(x=nkk, bins=best_bins, alpha=0.5, color = kallisto_color, label="Closest kallisto barcode")

    loc.set_xlabel('$\ell_1$ Distance', fontsize = fsize)
    loc.set_ylabel('Barcode counts', fontsize = fsize)
    loc.set_title('E', fontweight='bold', fontsize = fsize, loc = 'left' )
    loc.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = gridalpha)
    handles, labels = loc.get_legend_handles_labels()
    loc.legend(handles[::-1], labels[::-1])

    # F ###  kallisto t-SNE ##########
    loc = ax_f

    # plot a normal t-SNE for kallisto, business as usual
    loc.scatter(kallisto.obsm['TSNE10'][:,0], kallisto.obsm['TSNE10'][:,1],s =10, c = kallisto_color, alpha = 1, edgecolors = 'none', label = 'kallisto' )
    loc.set_ylabel(str(' t-SNE 2'), fontsize=fsize)
    loc.set_xlabel(str( 't-SNE 1'), fontsize=fsize)
    # loc.set_title('t-SNE on TSVD 10 components',fontweight='bold')
    loc.set_title('F', fontweight='bold', fontsize = fsize, loc = 'left' )
    loc.set_yticklabels([])
    loc.set_xticklabels([])
    loc.tick_params(axis=u'both', which=u'both',length=0)
    loc.legend()

    # G ### CellRanger t-SNE #############
    loc = ax_g

    # plot a normal t-SNE, business as usual
    sc = loc.scatter(tenx.obsm['TSNE10'][:,0], tenx.obsm['TSNE10'][:,1],s =10, c = tenx_color, alpha = 1, edgecolors = 'none', label = 'Cell Ranger' )
    loc.set_ylabel(str('t-SNE 2'), fontsize=fsize)
    loc.set_xlabel(str('t-SNE 1'), fontsize=fsize)
    loc.set_title('G', fontweight='bold', fontsize = fsize, loc = 'left' )
    loc.set_yticklabels([])
    loc.set_xticklabels([])
    loc.tick_params(axis=u'both', which=u'both',length=0)
    loc.legend()

    # H ##### DE PLOT #########
    loc = ax_h

    # load the csv file with the DE data, you'll need to change this path to match where it is
    df = pd.read_csv('/home/single_cell_analysis/brain_storm/gsea_bar/' + dataset_shortname + '.csv')

    # define some snazzy custom colors
    fold_change_change_colors = {
    '(4,5]':lighten_color(kallisto_color, 1.4),
    '(3,4]':lighten_color(kallisto_color, 1.1),
    '(2,3]':lighten_color(kallisto_color, 0.8),    
    '(1,2]':lighten_color(kallisto_color, 0.5),
    '(0,1]':lighten_color(kallisto_color, 0.2), 
    '(-1,0]':lighten_color(tenx_color, 0.2),
    '(-2,-1]':lighten_color(tenx_color, 0.5),
    '(-3,-2]':lighten_color(tenx_color, 0.8),
    '(-4,-3]':lighten_color(tenx_color, 1.1),
    '(-5,-4]':lighten_color(tenx_color, 1.4),
    }

    # For some datasets there are no DE genes, so we need this check to just write a text and make not plot
    if len(df)==0:
        loc.text(0.5*(1), 0.5*(1), 'No significant gene sets found ',
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=20, color='black',
            transform=loc.transAxes)
        loc.axis('off')

    # If there are DE genes, then we make a plot
    if len(df)>0:
        loc.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = 0.5)

        genesets = np.sort(df.gene_set.unique())
        changes = df.change.unique()
        changes =  list(reversed(changes))
        changes_dict = {}
        for change_interval in changes:
            changes_dict[change_interval] = []
            for geneset in genesets:
                ngenes_found = np.sum(df.loc[df['gene_set']==geneset]['change']==change_interval)
                changes_dict[change_interval].append(ngenes_found)
        ind = [x for x, _ in enumerate(genesets)]
        bottom_list = [0]*len(genesets)        
        for change_interval in changes_dict:
            plt.bar(ind, changes_dict[change_interval], bottom = bottom_list, label = change_interval, 
                    alpha = 0.5, width=0.8, color = fold_change_change_colors[change_interval])
            bottom_list = np.array(bottom_list) + np.array(changes_dict[change_interval])

        handles, labels = loc.get_legend_handles_labels()
        loc.legend(handles[::-1], labels[::-1], title='Fold change')

        loc.set_xticks(ticks=ind)
        loc.set_xticklabels( labels=genesets)
        loc.set_ylabel("Number of genes", fontsize=fsize)
        loc.set_xlabel("Gene set", fontsize=fsize)
        for label in loc.get_xmajorticklabels():
            label.set_rotation(30)
            label.set_horizontalalignment("right")

    loc.set_title('H', fontweight='bold', fontsize = fsize, loc = 'left' )

    # save the figure somewhere
    # plt.savefig(str('./figs/' + dataset_shortname + '.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(str('./figs/' + dataset_shortname + '.png'), dpi=300, bbox_inches='tight')

    print('-----------PLOT SAVED!---------')
    print(datetime.datetime.now())
    # plt.show()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="make plot")
    parser.add_argument("--ds", help="dataset name (i.e. SRR6998058_v2)")

    args = parser.parse_args()
    print("Loading files..")
    make_plot(str(args.ds))

    print("Done")
