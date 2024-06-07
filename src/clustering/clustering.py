import pandas

class Clustering:

    def __init__(self, matrix : pandas.DataFrame, k : int, max_iterations : int, centroid_labels : list, hvg):
        self.labels = self.kmeans_clustering(matrix, k, max_iterations, centroid_labels, hvg)
    
    def hvg(self, matrix : pandas.DataFrame, hvg : int):
        return matrix.std(0).sort_values(ascending = False).index.to_list()[0:hvg]

    def kmeans_clustering(self, matrix : pandas.DataFrame, k : int, max_iterations : int, centroid_labels : list, hvg) -> pandas.Series:
        if type(hvg) == list:
            hvg_matrix = matrix.loc[:, hvg]
        elif type(hvg) == int:
            hvg_matrix = matrix.loc[:, self.hvg(matrix, hvg)]
        elif type(hvg) == pandas.DataFrame:
            hvg_matrix = hvg
        centroids = hvg_matrix.loc[centroid_labels, :]
        centroids.index = [i for i in range(0, k)]

        for iter in range(0, max_iterations):

            distances = pandas.DataFrame()

            for i in range(0, len(centroids.index)):
                centroid = centroids.iloc[i, :]
                euclid_series = (hvg_matrix.sub(centroid, axis = 1) ** 2).sum(1) ** 1/2
                distances.insert(loc = i, column = i, value = euclid_series)

            labels = distances.idxmin(1)

            new_centroids = pandas.DataFrame(data = centroids, copy = True)
            for j in range(0, k):
                new_centroids.loc[j] = hvg_matrix.loc[labels == j].mean(0)

            if centroids.equals(new_centroids):
                break

            centroids = new_centroids
    
        return labels