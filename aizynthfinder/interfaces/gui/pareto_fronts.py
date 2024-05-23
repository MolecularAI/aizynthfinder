"""Module containing a GUI extension for plotting pareto fronts
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from IPython.display import display
from ipywidgets import Button, Checkbox, Dropdown, HBox, IntRangeSlider, Output, VBox
from aizynthfinder.analysis import (
    RouteCollection,
    RouteSelectionArguments,
    TreeAnalysis,
)
from aizynthfinder.interfaces.gui.utils import pareto_fronts_plot, route_display

if TYPE_CHECKING:
    from aizynthfinder.context.scoring.collection import ScorerCollection
    from aizynthfinder.interfaces.aizynthapp import AiZynthApp
    from aizynthfinder.search.andor_trees import AndOrSearchTreeBase
    from aizynthfinder.search.mcts import MctsSearchTree
    from aizynthfinder.utils.type_utils import Optional, StrDict, Union


class ParetoFrontsGUI:
    """GUI extension for plotting pareto fronts of routes.

    :param routes: the routes to plot
    :param scorers: the available scorers
    """

    def __init__(
        self,
        tree: Union[MctsSearchTree, AndOrSearchTreeBase],
        scorers: ScorerCollection,
    ) -> None:
        self._tree = tree
        self._scorers = scorers
        self._routes: Optional[RouteCollection] = None
        self._input: StrDict = dict()
        self._output: StrDict = dict()
        self._buttons: StrDict = dict()
        self._create_input()
        self._create_route_widgets()

    @classmethod
    def from_app(cls, app: AiZynthApp):
        """,
        Helper function to create a GUI from a GUI app interface

        :param app: the app to extract the routes from
        :return: the GUI object
        """
        if app.finder.tree is None:
            raise ValueError("Cannot initialize GUI from no tree")
        return ParetoFrontsGUI(app.finder.tree, app.finder.scorers)

    def _toggle_button(self, on_) -> None:
        for button in self._buttons.values():
            button.disabled = not on_

    def _create_input(self) -> None:
        self._input["scorer_1"] = Dropdown(
            options=self._scorers.names(),
            description="Scorer 1:",
            style={"description_width": "initial"},
        )
        self._input["scorer_2"] = Dropdown(
            options=self._scorers.names(),
            description="Scorer 2:",
            style={"description_width": "initial"},
        )
        self._input["select_routes"] = IntRangeSlider(
            value=[5, 25],
            min=1,
            max=50,
            step=1,
            description="Select #Routes:",
            disabled=False,
            continuous_update=True,
            orientation="horizontal",
            readout=True,
            readout_format="d",
            style={"description_width": "initial"},
        )
        self._input["select_all_routes"] = Checkbox(
            value=False,
            description="Select All Routes",
            style={"description_width": "initial"},
            layout={"justify": "left"},
        )
        selection_box = HBox(
            [
                self._input["select_routes"],
                self._input["select_all_routes"],
            ]
        )
        scorers_box = HBox(
            [
                self._input["scorer_1"],
                self._input["scorer_2"],
            ]
        )
        display(
            VBox([selection_box, scorers_box], layout={"border": "1px solid silver"})
        )

    def _create_route_widgets(self) -> None:
        self._buttons["show_routes"] = Button(description="Show Routes")
        self._buttons["show_routes"].on_click(self._on_display_routes_button_clicked)
        self._input["route"] = Dropdown(
            options=[],
            description="Routes: ",
        )
        self._input["route"].observe(self._on_change_route_option)
        display(
            HBox(
                [
                    self._buttons["show_routes"],
                    self._input["route"],
                ]
            )
        )

        self._output["pareto_fronts"] = Output(
            layout={"border": "1px solid silver", "width": "99%", "overflow": "auto"}
        )
        display(self._output["pareto_fronts"])

        self._output["routes"] = Output(
            layout={"border": "1px solid silver", "width": "99%", "overflow": "auto"}
        )
        display(self._output["routes"])

    def _on_display_routes_button_clicked(self, _) -> None:
        self._toggle_button(False)

        objectives = [self._input["scorer_1"].value, self._input["scorer_2"].value]
        mo_scorers = [self._scorers[name] for name in objectives]
        analysis = TreeAnalysis(self._tree, scorer=mo_scorers)

        config_selection = RouteSelectionArguments(
            nmin=self._input["select_routes"].value[0],
            nmax=self._input["select_routes"].value[1],
            return_all=self._input["select_all_routes"].value,
        )
        self._routes = RouteCollection.from_analysis(analysis, config_selection)
        self._routes.compute_scores(*self._scorers.objects())
        self._routes.make_images()

        self._input["route"].options = [
            f"Option {i}" for i, _ in enumerate(self._routes, 1)  # type: ignore
        ]

        with self._output["pareto_fronts"]:
            pareto_fronts_plot(self._routes)

        self._show_route(0)
        self._toggle_button(True)

    def _on_change_route_option(self, change) -> None:
        if change["name"] != "index":
            return
        self._show_route(self._input["route"].index)

    def _show_route(self, index) -> None:
        if self._routes is None:
            return
        route_display(index, self._routes, self._output["routes"])

    def _reset_widgets(self):
        self._output["pareto_fronts"].clear_output()
        self._output["routes"].clear_output()
        self._input["route"].options = []
