""" Module containing classes and routines for the GUI interface
"""
import tempfile
import signal
import subprocess
import os
import argparse
import logging

import ipywidgets as widgets
import jupytext
from ipywidgets import (
    HBox,
    Label,
    VBox,
    Text,
    Output,
    Checkbox,
    IntText,
    FloatText,
    Button,
    Dropdown,
    BoundedIntText,
    BoundedFloatText,
)
from rdkit import Chem
from IPython.display import display, HTML

from aizynthfinder.aizynthfinder import AiZynthFinder
from aizynthfinder.utils.logging import setup_logger


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
    :vartype finder: AiZynthFinder

    :param configfile: the path to yaml file with configuration
    :type configfile: str
    :param setup: if True will create and display the GUI on instatiation, defaults to True
    :type setup: bool, optional
    """

    def __init__(self, configfile, setup=True):
        setup_logger(logging.INFO)
        self.finder = AiZynthFinder(configfile=configfile)
        self._input = dict()
        self._output = dict()
        self._buttons = dict()
        if setup:
            self.setup()

    def setup(self):
        """
        Create the widgets and display the GUI.
        This is typically done on instatiation, but this method
        if for more advanced uses.
        """
        self._create_input_widgets()
        self._create_search_widgets()
        self._create_route_widgets()

    def _create_input_widgets(self):
        self._input["smiles"] = Text(description="SMILES", continuous_update=False)
        self._input["smiles"].observe(self._show_mol, names="value")
        display(self._input["smiles"])
        self._output["smiles"] = Output(
            layout={"border": "1px solid silver", "width": "50%", "height": "180px"}
        )
        display(self._output["smiles"])

        self._input["stocks"] = [
            Checkbox(value=True, description=key, layout={"justify": "left"})
            for key in self.finder.stock.available_stocks()
        ]
        box_stocks = VBox(
            [Label("Stocks")] + self._input["stocks"],
            layout={"border": "1px solid silver"},
        )

        self._input["policy"] = widgets.Dropdown(
            options=self.finder.policy.available_policies(),
            description="Neural Policy:",
            style={"description_width": "initial"},
        )

        max_time_box = self._make_slider_input("time_limit", "Time (min)", 1, 120)
        self._input["time_limit"].value = self.finder.config.time_limit / 60
        max_iter_box = self._make_slider_input(
            "iteration_limit", "Max Iterations", 100, 2000
        )
        self._input["iteration_limit"].value = self.finder.config.iteration_limit
        self._input["return_first"] = widgets.Checkbox(
            value=self.finder.config.return_first,
            description="Return first solved route",
        )
        vbox = VBox(
            [
                self._input["policy"],
                max_time_box,
                max_iter_box,
                self._input["return_first"],
            ]
        )
        box_options = HBox([box_stocks, vbox])

        self._input["C"] = FloatText(description="C", value=self.finder.config.C)
        self._input["max_transforms"] = BoundedIntText(
            description="Max steps for substrates",
            min=1,
            max=6,
            value=self.finder.config.max_transforms,
            style={"description_width": "initial"},
        )
        self._input["cutoff_cumulative"] = BoundedFloatText(
            description="Policy cutoff cumulative",
            min=0,
            max=1,
            value=self.finder.config.cutoff_cumulative,
            style={"description_width": "initial"},
        )
        self._input["cutoff_number"] = BoundedIntText(
            description="Policy cutoff number",
            min=1,
            max=1000,
            value=self.finder.config.cutoff_number,
            style={"description_width": "initial"},
        )
        self._input["exclude_target_from_stock"] = widgets.Checkbox(
            value=self.finder.config.exclude_target_from_stock,
            description="Exclude target from stock",
        )
        box_advanced = VBox(
            [
                self._input["C"],
                self._input["max_transforms"],
                self._input["cutoff_cumulative"],
                self._input["cutoff_number"],
                self._input["exclude_target_from_stock"],
            ]
        )

        children = [box_options, box_advanced]
        tab = widgets.Tab()
        tab.children = children
        tab.set_title(0, "Options")
        tab.set_title(1, "Advanced")
        display(tab)

    def _create_route_widgets(self):
        self._buttons["show_routes"] = Button(description="Show Reactions")
        self._buttons["show_routes"].on_click(self._on_display_button_clicked)
        self._input["route"] = Dropdown(options=[], description="Routes: ",)
        self._input["route"].observe(self._on_change_route_option)
        display(HBox([self._buttons["show_routes"], self._input["route"]]))

        self._output["routes"] = widgets.Output(
            layout={"border": "1px solid silver", "width": "99%"}
        )
        display(self._output["routes"])

    def _create_search_widgets(self):
        self._buttons["execute"] = Button(description="Run Search")
        self._buttons["execute"].on_click(self._on_exec_button_clicked)

        self._buttons["extend"] = widgets.Button(description="Extend Search")
        self._buttons["extend"].on_click(self._on_extend_button_clicked)
        display(HBox([self._buttons["execute"], self._buttons["extend"]]))

        self._output["tree_search"] = widgets.Output(
            layout={
                "border": "1px solid silver",
                "width": "99%",
                "height": "300px",
                "overflow": "auto",
            }
        )
        display(self._output["tree_search"])

    def _make_slider_input(self, label, description, min_val, max_val):
        label_widget = Label(description)
        slider = widgets.IntSlider(
            continuous_update=True, min=min_val, max=max_val, readout=False
        )
        self._input[label] = IntText(continuous_update=True, layout={"width": "80px"})
        widgets.link((self._input[label], "value"), (slider, "value"))
        return HBox([label_widget, slider, self._input[label]])

    def _on_change_route_option(self, change):
        if change["name"] != "index":
            return
        self._show_route(self._input["route"].index)

    def _on_exec_button_clicked(self, _):
        self._toggle_button(False)
        self._prepare_search()
        self._tree_search()
        self._toggle_button(True)

    def _on_extend_button_clicked(self, _):
        self._toggle_button(False)
        self._tree_search()
        self._toggle_button(True)

    def _on_display_button_clicked(self, _):
        self._toggle_button(False)
        self.finder.build_routes()
        self.finder.routes.make_images()
        self._input["route"].options = [
            f"Option {i}" for i, _ in enumerate(self.finder.routes, 1)
        ]
        self._show_route(0)
        self._toggle_button(True)

    def _prepare_search(self):
        self._output["tree_search"].clear_output()
        with self._output["tree_search"]:
            selected_stocks = [
                cb.description for cb in self._input["stocks"] if cb.value
            ]
            self.finder.stock.select_stocks(selected_stocks)
            self.finder.policy.select_policy(self._input["policy"].value)
            self.finder.config.update(
                **{
                    "C": self._input["C"].value,
                    "max_transforms": self._input["max_transforms"].value,
                    "cutoff_cumulative": self._input["cutoff_cumulative"].value,
                    "cutoff_number": int(self._input["cutoff_number"].value),
                    "return_first": self._input["return_first"].value,
                    "time_limit": self._input["time_limit"].value * 60,
                    "iteration_limit": self._input["iteration_limit"].value,
                    "exclude_target_from_stock": self._input[
                        "exclude_target_from_stock"
                    ].value,
                }
            )

            smiles = self._input["smiles"].value
            print("Setting target molecule with smiles: %s" % smiles)
            self.finder.target_smiles = smiles
            self.finder.prepare_tree()

    def _show_mol(self, change):
        self._output["smiles"].clear_output()
        with self._output["smiles"]:
            mol = Chem.MolFromSmiles(change["new"])
            display(mol)

    def _show_route(self, index):
        if index is None or index >= len(self.finder.routes):
            return

        node = self.finder.routes[index]["node"]
        state = node.state
        status = "Solved" if state.is_solved else "Not Solved"

        self._output["routes"].clear_output()
        with self._output["routes"]:
            display(HTML("<H2>%s" % status))
            display("Route Score: %0.3F" % state.score)
            display(HTML("<H2>Compounds to Procure"))
            display(state.to_image())
            display(HTML("<H2>Steps"))
            display(self.finder.routes[index]["image"])

    def _toggle_button(self, on):
        for button in self._buttons.values():
            button.disabled = not on

    def _tree_search(self):
        with self._output["tree_search"]:
            self.finder.tree_search(show_progress=True)
            print("Tree search completed.")


def _get_arguments():
    parser = argparse.ArgumentParser("aizynthapp")
    parser.add_argument(
        "--config", required=True, help="the filename of a configuration file"
    )
    parser.add_argument("--output", help="the filename of the Jupyter notebook")
    return parser.parse_args()


def main():
    """ Entry point for the aizynthapp command
    """
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
