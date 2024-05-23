""" Module containing classes and routines for the GUI interface
"""
import argparse
import logging
import os
import signal
import subprocess
import tempfile
from typing import TYPE_CHECKING

import ipywidgets as widgets
import jupytext
from IPython.display import HTML, display
from ipywidgets import (
    BoundedIntText,
    Button,
    Checkbox,
    Dropdown,
    HBox,
    IntText,
    Label,
    Output,
    SelectMultiple,
    Text,
    VBox,
)
from rdkit import Chem

from aizynthfinder.aizynthfinder import AiZynthFinder
from aizynthfinder.interfaces.gui.utils import pareto_fronts_plot, route_display
from aizynthfinder.utils.logging import setup_logger

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import StrDict


class AiZynthApp:
    """
    Interface class to be used in a Jupyter Notebook.
    Provides a basic GUI to setup and analyze the tree search.

    Should be instantiated with the path of a yaml file with configuration:

    .. code-block::

        from aizynthfinder.interfaces import AiZynthApp
        configfile = "/path/to/configfile.yaml"
        app = AiZynthApp(configfile)

    :ivar finder: the finder instance

    :param configfile: the path to yaml file with configuration
    :param setup: if True will create and display the GUI on instantiation, defaults to True
    """

    def __init__(self, configfile: str, setup: bool = True) -> None:
        setup_logger(logging.INFO)
        self.finder = AiZynthFinder(configfile=configfile)
        self._input: StrDict = dict()
        self._output: StrDict = dict()
        self._buttons: StrDict = dict()
        if setup:
            self.setup()

    def setup(self) -> None:
        """
        Create the widgets and display the GUI.
        This is typically done on instantiation, but this method
        is for more advanced uses.
        """
        self._create_input_widgets()
        self._create_search_widgets()
        self._create_route_widgets()

    def _create_input_widgets(self) -> None:
        self._input["smiles"] = Text(description="SMILES", continuous_update=False)
        self._input["smiles"].observe(self._show_mol, names="value")
        display(self._input["smiles"])
        self._output["smiles"] = Output(
            layout={"border": "1px solid silver", "width": "50%", "height": "180px"}
        )
        display(self._output["smiles"])

        self._input["stocks"] = [
            Checkbox(
                value=True,
                description=key,
                style={"description_width": "initial"},
                layout={"justify": "left"},
            )
            for key in self.finder.stock.items
        ]

        list_ = [Label("Limit atom occurrences")]
        self._input["stocks_atom_count_on"] = []
        self._input["stocks_atom_count"] = []
        current_criteria = self.finder.stock.stop_criteria.get("counts", {})
        for atom in ["C", "O", "N"]:
            chk_box = Checkbox(
                value=atom in current_criteria,
                description=atom,
                layout={"justify": "left", "width": "80px"},
                style={"description_width": "initial"},
            )
            self._input["stocks_atom_count_on"].append(chk_box)
            inpt = BoundedIntText(
                value=current_criteria.get(atom, 0),
                min=0,
                layout={"width": "80px"},
            )
            self._input["stocks_atom_count"].append(inpt)
            list_.append(HBox([chk_box, inpt]))
        box_stocks = VBox(
            [Label("Stocks")] + self._input["stocks"] + list_,
            layout={"border": "1px solid silver"},
        )

        init_value = (
            [self.finder.expansion_policy.items[0]]
            if self.finder.expansion_policy
            else []
        )
        self._input["policy"] = SelectMultiple(
            options=self.finder.expansion_policy.items,
            value=init_value,
            description="Expansion Policy:",
            style={"description_width": "initial"},
            rows=min(len(self.finder.expansion_policy.items) + 1, 4),
        )

        self._input["filter"] = SelectMultiple(
            options=self.finder.filter_policy.items,
            description="Filter Policy:",
            style={"description_width": "initial"},
            rows=min(len(self.finder.filter_policy.items) + 1, 4),
        )

        max_time_box = self._make_slider_input("time_limit", "Time (min)", 1, 120)
        self._input["time_limit"].value = self.finder.config.search.time_limit / 60
        max_iter_box = self._make_slider_input(
            "iteration_limit", "Max Iterations", 100, 2000
        )
        self._input["iteration_limit"].value = self.finder.config.search.iteration_limit
        self._input["return_first"] = widgets.Checkbox(
            value=self.finder.config.search.return_first,
            description="Return first solved route",
        )
        self._input["max_transforms"] = BoundedIntText(
            description="Max tree depth",
            min=1,
            max=20,
            value=self.finder.config.search.max_transforms,
            style={"description_width": "initial"},
        )

        rewards = self.finder.config.search.algorithm_config.get("search_rewards", [])
        objective1 = rewards[0] if len(rewards) > 0 else "state score"
        objective2 = rewards[1] if len(rewards) > 1 else "None"
        self._input["reward1"] = widgets.Dropdown(
            options=self.finder.scorers.names(),
            description="MCTS reward 1:",
            value=objective1,
            style={"description_width": "initial"},
        )

        self._input["reward2"] = widgets.Dropdown(
            options=["None"] + self.finder.scorers.names(),
            description="MCTS reward 2:",
            value=objective2,
            style={"description_width": "initial"},
        )

        vbox = VBox(
            [
                self._input["policy"],
                self._input["filter"],
                max_time_box,
                max_iter_box,
                self._input["return_first"],
                self._input["max_transforms"],
                self._input["reward1"],
                self._input["reward2"],
            ]
        )
        box_options = HBox([box_stocks, vbox])
        display(box_options)

    def _create_route_widgets(self) -> None:
        self._input["scorer"] = widgets.Dropdown(
            options=self.finder.scorers.names(),
            description="Reorder by:",
            style={"description_width": "initial"},
        )
        self._input["scorer"].disabled = True
        self._input["scorer"].observe(self._on_change_scorer)
        self._buttons["show_routes"] = Button(description="Show Routes")
        self._buttons["show_routes"].on_click(self._on_display_button_clicked)
        self._input["route"] = Dropdown(
            options=[],
            description="Routes: ",
        )
        self._input["route"].disabled = True
        self._input["route"].observe(self._on_change_route_option)
        display(
            HBox(
                [
                    self._buttons["show_routes"],
                    self._input["route"],
                    self._input["scorer"],
                ]
            )
        )

        self._output["pareto_fronts"] = widgets.Output(
            layout={"border": "1px solid silver", "width": "99%", "overflow": "auto"}
        )
        display(self._output["pareto_fronts"])
        self._output["routes"] = widgets.Output(
            layout={"border": "1px solid silver", "width": "99%"}
        )
        display(self._output["routes"])

    def _create_search_widgets(self) -> None:
        self._buttons["execute"] = Button(description="Run Search")
        self._buttons["execute"].on_click(self._on_exec_button_clicked)

        self._buttons["extend"] = widgets.Button(description="Extend Search")
        self._buttons["extend"].on_click(self._on_extend_button_clicked)
        display(HBox([self._buttons["execute"], self._buttons["extend"]]))

        self._output["tree_search"] = widgets.Output(
            layout={
                "border": "1px solid silver",
                "width": "99%",
                "height": "320px",
                "overflow": "auto",
            }
        )
        display(self._output["tree_search"])

    def _make_slider_input(self, label, description, min_val, max_val) -> HBox:
        label_widget = Label(description)
        slider = widgets.IntSlider(
            continuous_update=True, min=min_val, max=max_val, readout=False
        )
        self._input[label] = IntText(continuous_update=True, layout={"width": "80px"})
        widgets.link((self._input[label], "value"), (slider, "value"))
        return HBox([label_widget, slider, self._input[label]])

    def _on_change_route_option(self, change) -> None:
        if change["name"] != "index":
            return
        self._show_route(self._input["route"].index)

    def _on_change_scorer(self, change) -> None:
        if self.finder.routes is None or change["name"] != "index":
            return
        scorer = self.finder.scorers[self._input["scorer"].value]
        self.finder.routes.rescore(scorer)
        self._show_route(self._input["route"].index)

    def _on_exec_button_clicked(self, _) -> None:
        self._toggle_button(False)
        self._prepare_search()
        self._tree_search()
        self._toggle_button(True)

    def _on_extend_button_clicked(self, _) -> None:
        self._toggle_button(False)
        self._tree_search()
        self._toggle_button(True)

    def _on_display_button_clicked(self, _) -> None:
        self._toggle_button(False)
        rewards = self.finder.config.search.algorithm_config["search_rewards"]
        self.finder.build_routes(scorer=rewards)
        self.finder.routes.make_images()
        self.finder.routes.compute_scores(*self.finder.scorers.objects())
        self._input["route"].options = [
            f"Option {i}" for i, _ in enumerate(self.finder.routes, 1)  # type: ignore
        ]

        self._output["pareto_fronts"].clear_output()
        if len(rewards) > 1:
            with self._output["pareto_fronts"]:
                pareto_fronts_plot(self.finder.routes)

        self._show_route(0)
        self._input["scorer"].disabled = len(rewards) > 1
        self._input["route"].disabled = False
        self._toggle_button(True)

    def _prepare_search(self) -> None:
        self._output["tree_search"].clear_output()
        with self._output["tree_search"]:
            selected_stocks = [
                cb.description for cb in self._input["stocks"] if cb.value
            ]
            self.finder.stock.select(selected_stocks)
            atom_count_limits = {}
            for cb_input, value_input in zip(
                self._input["stocks_atom_count_on"], self._input["stocks_atom_count"]
            ):
                if cb_input.value:
                    atom_count_limits[cb_input.description] = value_input.value
            self.finder.stock.set_stop_criteria({"counts": atom_count_limits})
            self.finder.expansion_policy.select(self._input["policy"].value)
            if not self._input["filter"].value:
                self.finder.filter_policy.deselect()
            else:
                self.finder.filter_policy.select(self._input["filter"].value)

            self.finder.config.search.max_transforms = self._input[
                "max_transforms"
            ].value
            self.finder.config.search.return_first = self._input["return_first"].value
            self.finder.config.search.time_limit = self._input["time_limit"].value * 60
            self.finder.config.search.iteration_limit = self._input[
                "iteration_limit"
            ].value

            rewards = [self._input["reward1"].value]
            if self._input["reward2"].value != "None":
                rewards.append(self._input["reward2"].value)
            self.finder.config.search.algorithm_config["search_rewards"] = rewards
            if self.finder.config.search.algorithm != "mcts":
                print(
                    "Only MCTS algorithm is supported in GUI interface. Resetting to MCTS"
                )
                self.finder.config.search.algorithm = "mcts"

            smiles = self._input["smiles"].value
            print(f"Setting target molecule with smiles: {smiles}")
            self.finder.target_smiles = smiles
            self.finder.prepare_tree()

    def _show_mol(self, change) -> None:
        self._output["smiles"].clear_output()
        with self._output["smiles"]:
            mol = Chem.MolFromSmiles(change["new"])
            display(mol)

    def _show_route(self, index) -> None:
        route_display(index, self.finder.routes, self._output["routes"])

    def _toggle_button(self, on_) -> None:
        for button in self._buttons.values():
            button.disabled = not on_

    def _tree_search(self) -> None:
        with self._output["tree_search"]:
            self.finder.tree_search(show_progress=True)
            display(HTML("<b>Tree search completed!</b>"))


def _get_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser("aizynthapp")
    parser.add_argument(
        "--config", required=True, help="the filename of a configuration file"
    )
    parser.add_argument("--output", help="the filename of the Jupyter notebook")
    return parser.parse_args()


def main() -> None:
    """Entry point for the aizynthapp command"""
    args = _get_arguments()

    commands = "\n".join(
        [
            "from aizynthfinder.interfaces import AiZynthApp",
            f'configfile=r"{os.path.abspath(args.config)}"',
            "app = AiZynthApp(configfile)",
        ]
    )
    notebook = jupytext.reads(commands, fmt="py:percent")
    if args.output:
        filename = args.output
    else:
        _, filename = tempfile.mkstemp(suffix=".ipynb")
    jupytext.write(notebook, filename, fmt="ipynb")

    if args.output:
        print(f"Notebook saved to {filename}. It can now be open with Jupyter notebook")
        return

    try:
        proc = subprocess.Popen(f"jupyter notebook {filename}".split())
        proc.communicate()
    except KeyboardInterrupt:
        proc.send_signal(signal.SIGINT)


if __name__ == "__main__":
    main()
