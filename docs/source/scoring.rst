Scoring
=======

aizynthfinder is capable of scoring reaction routes, both in the form of ``Node`` objects when a search tree is available,
and in the form of ``ReactionTrees`` if post-processing is required.

Currently, there are a few scoring functions available

    * State score - a function of the number of precursors in stock and the length of the route
    * Number of reactions - the number of steps in the route
    * Number of pre-cursors - the number of pre-cursors in the route
    * Number of pre-cursors in stock - the number of the pre-cursors that are purchaseable


The *State score* is the score that is guiding the tree search in the :doc:`update phase <sequences>`, and 
this is not configurable. 

In the Jupyter notebook :doc:`GUI <gui>` one can choose to score the routes with any of the loaded the scorers. 

The three above scoring functions are loaded automatically when an ``aizynthfinder`` object is created.


Add new scoring functions
-------------------------


Additional scoring functions can be implemented by inheriting from the class ``Scorer`` in the ``aizynthfinder.context.scoring`` module.
The scoring class needs to implement the ``_score_node``, ``_score_reaction_tree`` and the ``__repr__`` methods.

This is an example of that.

.. code-block:: python

    from aizynthfinder.context.scoring import Scorer

    class DeltaNumberOfTransformsScorer(Scorer):

        def __repr__(self):
            return "delta number of transforms"

        def _score_node(self, node):
            return self._config.max_transforms - node.state.max_transforms

        def _score_reaction_tree(self, tree):
            return self._config.max_transforms - len(list(tree.reactions()))


This can then be added to the ``scorers`` attribute of an ``aizynthfinderfinder`` object. The ``scorers`` attribute is a collection
of ``Scorer`` objects.

For instance to use this in the Jupyter notebook GUI, one can do

.. code-block:: python

    from aizynthfinder.interfaces import AiZynthApp
    app = AiZynthApp("config_local.yml", setup=False)
    scorer = DeltaNumberOfTransformsScorer(app.finder.config)
    app.finder.scorers.load(scorer)
    app.setup()

