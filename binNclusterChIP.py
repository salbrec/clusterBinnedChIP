"""
Clustering ChIP-seq profiles and evaluating clusterin quality.
Profiles are provided in bed file-format, describing ChIP-seq peaks.
The distance used within the clustering is the Jaccard-Disance.
The distance is either computed based on the provided peak files using bedtools.
However, the script can also be used to bin the peak-profiles prior to clustering.
In our article, linked in the README, we describe that binning has a positive
impact on the clustering quality.

Run "python cluNbinChIP.py --help" to display a detailed help text on the console.
A comprehensive description is provided in the GitHub README that includes examples as well.

date:	January 2023
author:	Steffen Albrecht
"""

import pandas as pd
import numpy as np
import utils.utils as utils
from os.path import exists, abspath, realpath
from os import makedirs, getcwd
import matplotlib.pyplot as plt
from seaborn import clustermap
from sklearn.cluster import AgglomerativeClustering
from sklearn import metrics
import argparse
import pickle
from sys import argv
import time

class WrongUsageException(ValueError):
    pass

import warnings
warnings.filterwarnings("ignore")

out_dir = './example_CTCF/'

argsParser = argparse.ArgumentParser(description='Binning and Clustering for bed files - from ChIP-seq experiments')
argsParser.add_argument('--profiles', '-p', type=str, required=True, help='''Input table describing the bed-files for the
    analysis specified by ID, name, file path and label.''')
argsParser.add_argument('--type', '-t', type=str, default='peaks', choices=['peaks','bins'],
    help='''Profile type which is either "peaks" or "bins". When using bins, the
    parameters "--size" and "--genome" are required.''')
argsParser.add_argument('--size', '-s', type=str, default='None', help='''Must be used to specify the bin size
    when using "bins" as profile type. Sizes are given in kb or bp, for eexample "500bp" or  "2kb".''')
argsParser.add_argument('--genome', '-g', type=str, default='None', choices=['GRCh38','GRCm38'],
    help='''Must be used to specify the bin size. Sizes are given in kb or bp, for eexample "500bp" or  "2kb".''')
argsParser.add_argument('--output', '-o', type=str, default='./output/', help='''Specify the output folder.
    This folder folder will then contain the clustermaps, possibly also a subfolder for bins, etc.''')

args = argsParser.parse_args()
out_dir = args.output
type_descr = args.type
if not out_dir.endswith('/'):
    out_dir += '/'
profiles = pd.read_csv(args.profiles, sep='\t')
dist_tab_file = '%sdistances_on_%s.tsv'%(out_dir, args.type)
working_dir = getcwd() + '/'
script_dir = utils.get_path_to_filename(abspath(argv[0]))
messages = utils.get_messages()

print()
bins_map = {}
if args.type == 'bins':
    type_descr += '_%s_%s'%(args.genome, args.size)
    if args.size == 'None' or args.genome == 'None':
        message = 'The parameters "--size" and "--genome" must be used when using "bins" as profile-type!'
        raise WrongUsageException(message)
    bins_dir = out_dir + 'bins_%s_%s/'%(args.genome, args.size)
    if not exists(bins_dir):
        makedirs(bins_dir)
    dist_tab_file = '%sdistances_on_%s_%s_%s.tsv'%(out_dir, args.type, args.genome, args.size)
    if not exists(dist_tab_file):
        print(messages['comp_dists'])
        applied_binning = False
        for pi1, pr1 in profiles.iterrows():
            bins_file = '%s%s.pickle'%(bins_dir,pr1['ID'])
            if exists(bins_file):
                bins_map[pr1['ID']] = pickle.load(open(bins_file, 'rb'))
            else:
                if not applied_binning:
                    print(messages['comp_bins'])
                    applied_binning = True
                start_time = time.time()
                chr_sizes_file = '%sutils/%s.tsv'%(script_dir, args.genome)
                bins, max_bin_ID, position_map = utils.get_bins_from_file(pr1['filePath'],
                    chr_sizes_file, args.size)
                pickle.dump(bins, open(bins_file, 'wb'))
                bins_map[pr1['ID']] = bins
                print('Binning was applyed to %s. It took %.2f seconds.'%(pr1['filePath'],
                    time.time() - start_time))

distances_tab = {'pair':[], 'distance':[]}
print('\n\n')
if not exists(dist_tab_file):
    print(messages['comp_dists'])
    for pi1, pr1 in profiles.iterrows():
        for pi2, pr2 in profiles.iterrows():
            # the distances used are metrics and it we have that d(u,v) = d(v,u)
            # because of this, we can skip some distance computations
            if pi2 < pi1:
                continue
            pair = utils.get_ID_pair(pr1['ID'], pr2['ID'])
            d = None
            if args.type == 'peaks':
                d = utils.jaccard_dist_peaks(pr1['filePath'], pr2['filePath'])
            if args.type == 'bins':
                d = utils.jaccard_dist_sets(bins_map[pr1['ID']], bins_map[pr2['ID']])
            distances_tab['pair'].append( pair )
            distances_tab['distance'].append( d )
            if len(distances_tab['pair']) % 100 == 0:
                print('... %d distances were computed ...'%(len(distances_tab['pair'])))
    pd.DataFrame(distances_tab).to_csv(dist_tab_file, sep='\t', index=False)
    print('Done with computing distances, they are now available under:', dist_tab_file)
else:
    print('Distances have already been computed and are used from the file:', dist_tab_file)
    distances_tab = pd.read_csv(dist_tab_file, sep='\t')
distance_map = dict(zip( distances_tab['pair'], distances_tab['distance'] ))

# here the pre-computed distances are summarized in the right data structure
# as it is expected by the seaborn clustermap function
dists, indices = {}, []
for pi1, pr1 in profiles.iterrows():
    indices.append(pr1['name'])
    dists_temp_row = []
    for pi2, pr2 in profiles.iterrows():
        pair = utils.get_ID_pair(pr1['ID'], pr2['ID'])
        dists_temp_row.append( distance_map[pair] )
    dists[pr1['name']] = dists_temp_row
dists = pd.DataFrame(dists, index=indices)

# creating the clustermap based on single-linkage agglomerative clustering
plt.figure(figsize=(4, 4))
cm = clustermap(dists, linewidth=3, method='single',
    cmap='Blues_r', col_cluster=True, row_cluster=True,
    xticklabels=True, yticklabels=True)
# a png is created for the clustermap
cluster_map_file = '%sclustermap_%s.png'%(out_dir,type_descr)
plt.savefig(cluster_map_file)
print('\nClustermap was created and saved into', cluster_map_file)

# depending on the availability of labels for the given profiles,
# the clustering is done again but now validated with common evaluation metrics
print()
eval_clustering = False
if not 'label' in profiles.columns:
    print('Skipping the clustering evaluation because profiles are not labeled.')
    print('(the provided table in "%s" has no column "label")\n\n'%(args.profiles))
elif len(set(profiles['label'])) > 1:
    eval_clustering = True
else:
    print('Skipping the clustering evaluation because all profiles are from the same label.\n\n')

if eval_clustering:
    label_names = list(profiles['label'])
    uniq_labels = list(set(label_names))
    true_labels = [ uniq_labels.index(ln) for ln in label_names ]

    agg = AgglomerativeClustering(n_clusters=len(uniq_labels),
        affinity='precomputed', linkage='single')
    pred_labels = agg.fit_predict(np.array(dists))

    ami = metrics.adjusted_mutual_info_score(true_labels, pred_labels)
    ari = metrics.adjusted_rand_score(true_labels, pred_labels)
    hos = metrics.homogeneity_score(true_labels, pred_labels)
    cos = metrics.completeness_score(true_labels, pred_labels)

    print('''Applied agglomerative single-linkage clustering and evaluated its quality.
The following measures represent the clustering quality:\n''')
    print('\tAdjusted Mutual Information: %.5f'%(ami))
    print('\tAdjusted Rand Index:         %.5f'%(ari))
    print('\tHomogeneity Score:           %.5f'%(hos))
    print('\tCompleteness Score:          %.5f'%(cos))
    print('\n')









