from aizynthfinder.analysis.routes import RouteCollection
from aizynthfinder.interfaces.gui.utils import pareto_fronts_plot, route_display
from aizynthfinder.interfaces.gui.pareto_fronts import ParetoFrontsGUI


def test_plot_pareto_fronts(setup_mo_scorer, setup_analysis, default_config):
    analysis, _ = setup_analysis(scorer=setup_mo_scorer(default_config))
    routes = RouteCollection.from_analysis(analysis)

    pareto_fronts_plot(routes)


def test_display_route(setup_analysis, mocker):
    display_patch = mocker.patch("aizynthfinder.interfaces.gui.utils.display")
    mocked_widget = mocker.MagicMock()
    analysis, _ = setup_analysis()
    routes = RouteCollection.from_analysis(analysis)
    routes.make_images()

    route_display(0, routes, mocked_widget)
    display_patch.assert_called()


def test_pareto_gui(setup_analysis, mocker, default_config):
    display_patch = mocker.patch("aizynthfinder.interfaces.gui.pareto_fronts.display")
    analysis, _ = setup_analysis()
    routes = RouteCollection.from_analysis(analysis)

    ParetoFrontsGUI(routes.reaction_trees, default_config.scorers)

    display_patch.assert_called()
