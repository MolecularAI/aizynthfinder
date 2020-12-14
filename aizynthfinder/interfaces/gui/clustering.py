""" Module containing a GUI extension for clustering
"""
import matplotlib.pylab as plt
import numpy as np
from ipywidgets import Output, Label, BoundedIntText, Button, HBox, Tab
from IPython.display import display
from scipy.cluster.hierarchy import dendrogram

from aizynthfinder.utils.route_clustering import ClusteringHelper


class ClusteringGui:
    """
    GUI extention to cluster routes

    :param routes: the routes to cluster
    :type routes: RouteCollection
    :param content: what to cluster on
    :type content: TreeContent, optional
    """

    def __init__(self, routes, content="both"):
        self._routes = routes
        self._routes.distance_matrix(content=content, recreate=True)
        self._input = dict()
        self._output = dict()
        self._buttons = dict()
        self._create_dendogram()
        self._create_input()
        self._create_output()

    @classmethod
    def from_app(self, app, content="both"):
        """
        Helper function to create a GUI from a GUI app interface

        :param app: the app to extract the routes from
        :type app: AiZynthApp
        :param content: what to cluster on
        :type content: TreeContent, optional
        :return: the GUI object
        :rtype: ClusteringGUI
        """
        return ClusteringGui(app.finder.routes, content)

    def _create_dendogram(self):
        dend_out = Output(
            layout={"width": "99%", "height": "310px", "overflow_y": "auto"}
        )
        with dend_out:
            print("This is the hierarchy of the routes")
            fig = plt.Figure()
            dendrogram(
                ClusteringHelper(self._routes.distance_matrix()).linkage_matrix(),
                color_threshold=0.0,
                labels=np.arange(1, len(self._routes) + 1),
                ax=fig.gca(),
            )
            fig.gca().set_xlabel("Route")
            fig.gca().set_ylabel("Distance")
            display(fig)
        display(dend_out)

    def _create_input(self):
        self._input["number_clusters"] = BoundedIntText(
            continuous_update=True,
            min=1,
            max=len(self._routes) - 1,
            layout={"width": "80px"},
        )
        self._buttons["cluster"] = Button(description="Cluster")
        self._buttons["cluster"].on_click(self._on_cluster_button_clicked)
        box = HBox(
            [
                Label("Number of clusters to make"),
                self._input["number_clusters"],
                self._buttons["cluster"],
            ]
        )
        display(box)
        help_out = Output()
        with help_out:
            print(
                "Optimization is carried out if the number of given clusters are less than 2"
            )
        display(help_out)

    def _create_output(self):
        self._output["clusters"] = Tab()
        display(self._output["clusters"])

    def _on_cluster_button_clicked(self, _):
        self._buttons["cluster"].enabled = False
        self._routes.cluster(self._input["number_clusters"].value)
        self._buttons["cluster"].enabled = True

        outputs = []
        for i, cluster in enumerate(self._routes.clusters):
            output = Output(
                layout={
                    "border": "1px solid silver",
                    "width": "99%",
                    "height": "500px",
                    "overflow_y": "auto",
                }
            )
            with output:
                for image in cluster.images:
                    print(f"Route {self._routes.images.index(image)+1}")
                    display(image)
            outputs.append(output)
            self._output["clusters"].set_title(i, f"Cluster {i+1}")
        self._output["clusters"].children = outputs
