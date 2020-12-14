""" Module containing a class to help out with clustering
"""
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score


class ClusteringHelper:
    """
    A helper class to perform clustering of items
    based on a pre-computed distance matrix.

    The clustering is the agglomerative clustering
    algorithm implemented in Scikit-learn.
    """

    def __init__(self, distances):
        self._distances = distances
        self._model = None
        self.optimization_scores = []

    @property
    def labels(self):
        return self._model.labels_ if self._model else None

    def fixed_clustering(self, n_clusters, **kwargs):
        """
        Make a fixed number of clusters.

        Additional arguments to the clustering algorithm
        can be passed in as key-word arguments

        :param n_clusters: the desired number of clusters
        :type n_clusters: int
        :return: the cluster index for each observation
        :rtype: numpy.ndarray
        """
        if "linkage" not in kwargs:
            kwargs["linkage"] = "single"
        kwargs["affinity"] = "precomputed"
        kwargs["n_clusters"] = n_clusters
        model = AgglomerativeClustering(**kwargs)
        model.fit(self._distances)
        self._model = model
        return model.labels_

    def linkage_matrix(self, **kwargs):
        """
        Compute the linkage matrix.

        Additional arguments to the clustering algorithm
        can be passed in as key-word arguments

        :return: the linkage matrix
        :rtype: numpy.ndarray
        """
        if "linkage" not in kwargs:
            kwargs["linkage"] = "single"
        kwargs["affinity"] = "precomputed"
        kwargs["n_clusters"] = None
        kwargs["distance_threshold"] = 0.0
        model = AgglomerativeClustering(**kwargs)
        model.fit(self._distances)

        self._model = model
        counts = np.zeros(len(model.distances_))
        matrix = np.column_stack([model.children_, model.distances_, counts])
        return matrix.astype(float)

    def optimize(self, max_clusters=5, **kwargs):
        """
        Optimize the number of cluster based  Silhouette metric.

        Additional arguments to the clustering algorithm
        can be passed in as key-word arguments

        :param max_clusters: the maximum number of clusters to consider
        :type max_clusters: int, optional
        :return: the cluster index for each observation
        :rtype: numpy.ndarray
        """

        max_score = None
        best_size = None
        self.optimization_scores = []
        for n_clusters in range(2, min(max_clusters + 1, len(self._distances))):
            clusters = self.fixed_clustering(n_clusters, **kwargs)
            score = silhouette_score(self._distances, clusters, metric="precomputed")
            self.optimization_scores.append(score)
            if best_size is None or score > max_score:
                max_score = score
                best_size = n_clusters

        return self.fixed_clustering(best_size, **kwargs)

    @staticmethod
    def cluster(distances, n_clusters, **kwargs):
        """
        Cluster items based on a pre-computed distance matrix using a hierarchical clustering.

        :param distances: the distance matrix
        :type distances: numpy.ndarray
        :param n_clusters: the desired number of clusters
        :type n_clusters: int
        :return: the cluster index for each observation
        :rtype: numpy.ndarray
        """
        helper = ClusteringHelper(distances)

        max_clusters = kwargs.pop("max_clusters", 5)

        if n_clusters < 2:
            return helper.optimize(max_clusters=max_clusters, **kwargs)

        return helper.fixed_clustering(n_clusters, **kwargs)
