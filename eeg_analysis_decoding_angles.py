#!/usr/bin/env python
# coding: utf-8

# In[1]:


# this notebook is for decoding the orientation of difference between left and right stimulus


# In[2]:


# load necessary libraries for decoding(categorical)
import numpy as np
import pandas as pd
import pickle
import os.path
import matplotlib.pyplot as plt
from joblib import Memory
from joblib import Parallel, delayed
import joblib
from scipy import stats
from pathlib import Path
from sklearn.pipeline import Pipeline
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import ShuffleSplit, cross_val_score
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

# import mne
import mne
from mne.datasets import sample
from mne.decoding import (SlidingEstimator, GeneralizingEstimator, Scaler,
                          cross_val_multiscore, LinearModel, get_coef,
                          Vectorizer, CSP)
from mne import (read_epochs, combine_evoked, write_evokeds)
from mne.stats import bonferroni_correction, fdr_correction
from mne import io


# In[3]:


# load the basic folder names and the location of the eeg datasets
folder_name_eeg = "E:\\SUMO_further_data_pack_zx\\N2pc_IEM\\new_results\\eeg_before_IEM\\"
folder_name_beh = "E:\\SUMO_further_data_pack_zx\\N2pc_IEM\\new_results\\eeg_beh\\"
subject_name = ['SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',
                'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120', 
                'SUMO_3001', 'SUMO_3015', 'SUMO_3017']
individual_name_eeg_tail = '_before_iem_stim.set'
individual_name_beh_tail = '_beh.csv'


# In[ ]:


# preallocate the space for group data
score_mean_group = np.empty([len(subject_name), 1000])
    
# load csv file and the eeg data
for s in subject_name:
    eeg_file_name = str(folder_name_eeg + s + individual_name_eeg_tail)
    beh_file_name = str(folder_name_beh + s + individual_name_beh_tail)

    # load eeg file to the mne workspace
    eeg = mne.read_epochs_eeglab(eeg_file_name)

    # load the csv beh file to the df
    df = pd.read_csv(beh_file_name)

    # copy the df to df2 for modification
    df2 = df.copy()
    df2 = df2[['type', 'epoch', 'orientdiff']]
    df2 = df2[df2['type'] == 'S 48'].reset_index()
    
    # do the k-means clustering to find the class numbers
    X = np.array(df2.orientdiff).reshape(-1, 1)
    # by using the sse and the elbow method, I pick the n_cluster equals to 5(see supplementaries)
    n_clusters = 5
    kmeans = KMeans(n_clusters=n_clusters, n_init=100, max_iter=3000).fit(X)
    label = kmeans.labels_
    
    # change from 1-len(label) to 0 - len(label)-1(from MATLAB to python)
    df2.epoch = df.epoch - 1
    df2['label'] = label
    
    # load parameters for later decoding
    number_of_cpu = joblib.cpu_count()
    clf = make_pipeline(Scaler(eeg.info),
                        Vectorizer(),
                        LogisticRegression(solver='lbfgs'))
    # CSP
    csp = CSP(n_components=5, norm_trace=False, cov_est='concat', cov_method_params = dict)
    clf_csp = make_pipeline(csp, LinearModel(LogisticRegression(solver='lbfgs')))
    # Time_decoding
    clf_1 = make_pipeline(StandardScaler(), LogisticRegression(solver='lbfgs'))
    time_decode = SlidingEstimator(clf_1, n_jobs=number_of_cpu, scoring='roc_auc', verbose=True)

    # pre-allocate the space
    score_mean = np.empty([n_clusters, len(eeg.times)])
    #print(df2.head())
    for i in list(range(n_clusters)):
        df3 = df2.copy()
        df3.loc[df2['label'] == i, 'label'] = 1
        df3.loc[df2['label'] != i, 'label'] = 0
        #print(df3.head())
        eeg2 = eeg.get_data()
        score = cross_val_multiscore (time_decode, eeg2, df3.label, cv = 10, n_jobs=number_of_cpu)
        score_mean[i] = np.mean(score, axis=0)

    p = np.mean(score_mean, axis=0)
    
    score_mean_group[subject_name.index(s)] = p


# In[ ]:


q = np.mean(score_mean_group, axis=0)
plt.figure(figsize=(10, 6), dpi=80)
plt.plot(eeg.times, q, 'b', label='Scene_post')
#plt.errorbar(eeg.times, q_fil, er)
xmin, xmax = plt.xlim()
# plt.hlines(0.5, xmin, xmax, linestyle='--', colors='k',
#            label='chance', linewidth=2)
#plt.hlines(threshold_fdr_f_p, xmin, xmax, linestyle='--', colors='b',
           #label='p=0.05 (FDR)', linewidth=2)
plt.legend()
plt.xlabel("Time (sec)")
plt.ylabel("AUC")
plt.savefig("stimulus_angles.pdf", format="pdf", bbox_inches="tight")
plt.savefig("stimulus_angles.png", format="png", bbox_inches="tight")
plt.show()

np.savetxt("stimulus_angles_group_score.csv", score_mean_group, delimiter=",")



# In[ ]:


# # Supplementary: show the SSE to see what the number of cluster can be most recommended
# import collections

# for s in subject_name:
#     beh_file_name = str(folder_name_beh + s + individual_name_beh_tail)
#     # load the csv beh file to the df
#     df = pd.read_csv(beh_file_name)

#     # copy the df to df2 for modification
#     df2 = df.copy()
#     df2 = df2[['type', 'epoch', 'orientdiff']]
#     df2 = df2[df2['type'] == 'S 48'].reset_index()
    
#     X = np.array(df2.orientdiff).reshape(-1, 1)
#     n_clusters = np.array(range(1,20))
#     sse = []
#     for n in n_clusters:
#         kmeans = KMeans(n_clusters=n, n_init=100, max_iter=3000).fit(X)
#         sse.append(kmeans.inertia_)
#         counter=collections.Counter(kmeans.labels_)
#         # show the number of trials in each cluster
#         print(n, counter)

#     plt.plot(range(1, 20), sse)
#     plt.xticks(range(1, 20))
#     plt.title(s)
#     plt.xlabel("Number of Clusters")
#     plt.ylabel("SSE")
#     plt.show()

