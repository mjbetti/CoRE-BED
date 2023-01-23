from dask import dataframe as dd
import dask
import dask.array as da
import dask_ml.datasets
import dask_ml.cluster
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar
ProgressBar().register()
import pickle
from kmodes.kmodes import KModes
import pandas as pd
import numpy as np
from sklearn import tree
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.tree import export_graphviz
from six import StringIO
from IPython.display import Image
import pydotplus
from sklearn.tree import export_text
from sklearn.decomposition import PCA
from sklearn.tree import export_text

#Open the expression data as a pandas data frame
path = "/home/bettimj/gamazon_rotation/mod_core-bed/epigenomes_expanded/kmeans_training_and_test_epigenome.txt"

print("Importing input file...")
file = dd.read_csv(path, sep = "\t")
#file = file.sample(frac=1)
filt = file.drop(['chrom', 'start', 'end', 'tss_5k1k', 'tissue'], axis=1)
pandas_df = filt.compute()
pandas_df_uniq = pandas_df.drop_duplicates()
print("Converting to array...")
#array = filt.to_dask_array(lengths=True)
np_df = pandas_df_uniq.to_numpy()

#wcss = []

print("Processing " + str(9) + " clusters...")
clustering = KModes(n_clusters = 9, init = "Cao", verbose = 1, n_init = 10)
clusters = clustering.fit_predict(np_df)

#PCA
##3D plot of PC1, PC2, and PC3
pca = PCA(3)
plot_columns = pca.fit_transform(pandas_df_uniq)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
for label in np.unique(clusters):
	ax.scatter(plot_columns[clusters == label, 0], plot_columns[clusters == label, 1], plot_columns[clusters == label, 2], label = label)
ax.legend()
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
fig.savefig("core_bed_kmodes_pca_k9.png", dpi = 300)

##2D plot of PC1 and PC2
pca = PCA(3)
plot_columns = pca.fit_transform(pandas_df_uniq)

fig = plt.figure()
ax = fig.add_subplot()
for label in np.unique(clusters):
	ax.scatter(plot_columns[clusters == label, 0], plot_columns[clusters == label, 1], label = label)
ax.legend()
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
fig.savefig("core_bed_kmodes_pca_k9_pc1_pc2.png", dpi = 300)

##2D plot of PC1 and PC3
fig = plt.figure()
ax = fig.add_subplot()
for label in np.unique(clusters):
	ax.scatter(plot_columns[clusters == label, 0], plot_columns[clusters == label, 2], label = label)
ax.legend()
ax.set_xlabel('PC1')
ax.set_ylabel('PC3')
fig.savefig("core_bed_kmodes_pca_k9_pc1_pc3.png", dpi = 300)

##2D plot of PC2 and PC3
fig = plt.figure()
ax = fig.add_subplot()
for label in np.unique(clusters):
	ax.scatter(plot_columns[clusters == label, 1], plot_columns[clusters == label, 2], label = label)
ax.legend()
ax.set_xlabel('PC2')
ax.set_ylabel('PC3')
fig.savefig("core_bed_kmodes_pca_k9_pc2_pc3.png", dpi = 300)

#Decision tree training
labels_df = clusters
labels_df_comp = labels_df.tolist()

filt_comp = pandas_df_uniq

x = filt_comp.to_numpy().tolist()
y = labels_df_comp

X_train, test_x, y_train, test_lab = train_test_split(x, y, test_size = 0.1, random_state = 42)

depths = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
feature_cols = ["tss_25k25k", "27ac", "27me3", "4me1", "4me3", "atac", "ctcf", "dnase", "ep300", "h2afz", "h3k36me3", "h3k4me2", "h3k79me2", "h3k9ac", "h3k9me3", "h4k20me1", "polr2a", "rad21", "smc3"]

for depth in depths:
	print(depth)
	clf = tree.DecisionTreeClassifier(random_state = 42, criterion = "entropy", max_depth = depth)
	clf = clf.fit(X_train, y_train)

	y_pred = clf.predict(test_x)
	print("Accuracy:",metrics.accuracy_score(test_lab, y_pred))
	acc = "Accuracy:" + str(metrics.accuracy_score(test_lab, y_pred))
	f = open("core_bed_tree_depth" + str(depth) + ".txt", "w")
	f.write(acc)

	pickle.dump(clf, open("decision_tree_k7_entropy.depth" + str(depth) + ".pickle", 'wb'))

	r = export_text(clf, feature_names=["tss_25k25k", "27ac", "27me3", "4me1", "4me3", "atac", "ctcf", "dnase", "ep300", "h2afz", "h3k36me3", "h3k4me2", "h3k79me2", "h3k9ac", "h3k9me3", "h4k20me1", "polr2a", "rad21", "smc3"], max_depth = depth)
	print(r)
	f = open("core_bed_tree_depth" + str(depth) + ".txt", "a")
	f.write(r)
	f.close()