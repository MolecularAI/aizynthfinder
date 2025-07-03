"""Module containing utility functions for GUI.
"""

from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING

import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from IPython.display import HTML, display
from paretoset import paretorank

if TYPE_CHECKING:
    from aizynthfinder.analysis.routes import RouteCollection
    from aizynthfinder.utils.type_utils import List, Optional


def pareto_fronts_plot(
    routes: RouteCollection,
) -> None:
    """Plot the pareto front(s).

    :param routes: the route collection to plot as Pareto fronts
    """
    scorer_names = list(routes.scores[0].keys())
    scores = np.array(
        [[score_dict[name] for name in scorer_names] for score_dict in routes.scores]
    )
    direction_arr = np.repeat("max", len(scorer_names))
    pareto_ranks = paretorank(scores, sense=direction_arr, distinct=False)

    pareto_fronts = pd.DataFrame(scores, columns=scorer_names)
    pareto_fronts.loc[:, "pareto_rank"] = pareto_ranks
    pareto_fronts.loc[:, "route"] = np.arange(1, scores.shape[0] + 1)
    pareto_fronts_unique = pareto_fronts.drop_duplicates(scorer_names)
    # Apply the default theme
    sns.set_theme()
    fig = sns.relplot(
        data=pareto_fronts_unique,
        x=scorer_names[0],
        y=scorer_names[1],
        hue="pareto_rank",
        kind="line",
        markers=True,
    )
    fig.set(
        xlabel=f"{scorer_names[0]}",
        ylabel=f"{scorer_names[1]}",
        title="Pareto Fronts",
    )

    objectives2solutions = defaultdict(list)
    for _, row in pareto_fronts.iterrows():
        x_val = row[scorer_names[0]]
        y_val = row[scorer_names[1]]
        objectives2solutions[(x_val, y_val)].append(int(row["route"]))

    # Add route label on pareto line(s)
    for _, row in pareto_fronts_unique.iterrows():
        x_val = row[scorer_names[0]]
        y_val = row[scorer_names[1]]
        point_val = f"Option {_values_to_string(objectives2solutions[(x_val, y_val)])}"
        plt.text(x=x_val, y=y_val, s=point_val, size=10)
    plt.show()


def route_display(
    index: Optional[int], routes: RouteCollection, output_widget: widgets.Output
) -> None:
    """
    Display a route with auxillary information in a widget

    :param index: the index of the route to display
    :param routes: the route collection
    :param output_widget: the widget to display the route on
    """
    if index is None or routes is None or index >= len(routes):
        return

    route = routes[index]
    state = route["node"].state
    status = "Solved" if state.is_solved else "Not Solved"

    output_widget.clear_output()
    with output_widget:
        display(HTML(f"<H2>{status}"))
        table_content = "".join(
            (
                f"<tr><td>{name}</td><td>{score:.4f}</td></tr>"
                if isinstance(score, (float, int))
                else f"<tr><td>{name}</td><td>({', '.join(f'{val:.4f}' for val in score)})</td></tr>"
            )
            for name, score in route["all_score"].items()
        )
        display(HTML(f"<table>{table_content}</table>"))
        display(HTML("<H2>Compounds to Procure"))
        display(state.to_image())
        display(HTML("<H2>Steps"))
        display(routes[index]["image"])


def _values_to_string(vals: List[int]) -> str:
    """
    Given a list of integers, produce a nice-looking string

    Example:
         >> _values_to_string([1,2,3,10,11,19])
         "1-3, 10-11, 19"
    """
    groups = []
    start = end = vals[0]
    for prev_val, val in zip(vals[:-1], vals[1:]):
        if prev_val == val - 1:
            end += 1
        else:
            groups.append((start, end))
            start = end = val
    groups.append((start, end))

    group_strs = []
    for start, end in groups:
        if start == end:
            group_strs.append(str(start))
        else:
            group_strs.append(f"{start}-{end}")
        if len(group_strs) % 5 == 0:
            group_strs[-1] = "\n" + group_strs[-1]
    return ", ".join(group_strs)
