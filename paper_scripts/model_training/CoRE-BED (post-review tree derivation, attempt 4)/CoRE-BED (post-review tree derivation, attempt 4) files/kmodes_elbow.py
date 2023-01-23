from dask import dataframe as dd
import dask.array as da
import dask_ml.datasets
import dask_ml.cluster
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar
ProgressBar().register()
import pickle
from kmodes.kmodes import KModes
import pandas as pd
import seaborn as sns

#Open the expression data as a pandas data frame
path = "/home/bettimj/gamazon_rotation/mod_core-bed/epigenomes_expanded/kmeans_elbow_epigenome.txt"

print("Importing input file...")
file = dd.read_csv(path, sep = "\t")
#file = file.sample(frac=1)
filt = file.drop(['chrom', 'start', 'end', 'tss_5k1k', 'tissue'], axis=1)
pandas_df = filt.compute()
pandas_df_uniq = pandas_df.drop_duplicates()
print("Converting to array...")
#array = filt.to_dask_array(lengths=True)
np_df = pandas_df_uniq.to_numpy()

wcss = []
ks = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

for i in ks:
	print("Processing " + str(i) + " clusters...")
	clustering = KModes(n_clusters = i, init = "Cao", verbose = 1, n_init = 10)
	clustering.fit_predict(np_df)
	wcss.append(clustering.cost_)
        
#Plot the data
plt = sns.lineplot(x = ks, y = wcss)
plt.set_xticks(range(1,21)) # <--- set the ticks first
plt.set_xticklabels(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'])
fig = plt.get_figure()
fig.savefig("core_bed_kmodes_elbow_k1_20.png", dpi = 300)