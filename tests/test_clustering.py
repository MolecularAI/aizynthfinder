import numpy as np
import pytest

from aizynthfinder.utils.route_clustering.clustering import ClusteringHelper


distance_matrix = np.array(
    [
        [0.0, 1.0, 20.0, 15.0],
        [1.0, 0.0, 17.0, 12.0],
        [20.0, 17.0, 0.0, 5.0],
        [15.0, 12.0, 5.0, 0.0],
    ]
)


def test_create_clustering_helper():
    helper = ClusteringHelper(distance_matrix)

    with pytest.raises(ValueError):
        _ = helper.labels


def test_make_two_clusters():
    helper = ClusteringHelper(distance_matrix)

    labels = helper.fixed_clustering(2)

    assert labels[0] == labels[1]
    assert labels[0] != labels[2]
    assert labels[2] == labels[3]


def test_optimize_clusters():
    helper = ClusteringHelper(distance_matrix)

    labels = helper.optimize()

    assert max(labels) == 1
    assert labels[0] == labels[1]
    assert labels[0] != labels[2]
    assert labels[2] == labels[3]
    assert len(helper.optimization_scores) == 2


def test_clustering_helper():
    labels1 = ClusteringHelper.cluster(distance_matrix, 2)
    labels2 = ClusteringHelper.cluster(distance_matrix, 1)

    assert list(labels1) == list(labels2)


def test_linkage_matrix():
    helper = ClusteringHelper(distance_matrix)

    matrix = helper.linkage_matrix()

    assert len(matrix) == 3
