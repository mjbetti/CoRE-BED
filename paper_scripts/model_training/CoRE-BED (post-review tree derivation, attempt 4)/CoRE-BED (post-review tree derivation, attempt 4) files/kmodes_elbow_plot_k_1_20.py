import seaborn as sns

wcss = [101391.0, 87976.0, 81553.0, 78041.0, 76294.0, 75216.0, 73862.0, 73927.0, 70504.0, 69372.0, 67488.0, 66438.0, 66046.0, 65067.0, 64851.0, 64270.0, 63313.0, 62252.0, 61335.0, 60892.0]

ks = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
plt = sns.lineplot(x = ks, y = wcss)
plt.set_xticks(range(1,21)) # <--- set the ticks first
plt.set_xticklabels(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'])
plt.set_xlabel('Number of clusters')
plt.set_ylabel('Cost')
fig = plt.get_figure()
fig.savefig("core_bed_kmodes_elbow_k1_20.png", dpi = 300)
