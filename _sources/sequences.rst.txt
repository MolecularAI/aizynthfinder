Sequences
=========

This page shows some sequence diagrams to aid in the understanding of how information
is passed between different objects in the Monte Carlo tree search.
The sequences are simplified, but explains the overall picture.
The flow of information / method call should be read top-down.

Analysis / post-processing
--------------------------

This sequence explains how the ``AiZynthFinder`` object exports the top-ranked reaction tree 
as a JSON. Note, this is only one possible sequence for analysis of the trees.

.. image:: analysis-seq.png


Monte Carlo tree search
-----------------------

This sequence explains how the Monte Carlo tree search is carried out by the ``AiZynthFinder`` object.

.. image:: treesearch-seq.png

The following text explains what is executed at each iteration of the tree search (the outer loop in the ``one_iteration()`` method of the ``MctsSearchTree`` class).

First, a leaf is selected using the ``select_leaf()`` method of the ``MctsSearchTree`` class. This is called the `Selection` phase in the literature and will pick the most promising leaf to continue the search from. In the first iteration, this is simply the root node. For the rest of the iterations, the algorithm will execute the following:

1.	Set the current node to the root
2.	Loop while the current is expanded, and the state of the current node is not solved

  a. Select the most promising child of the current node by calling the ``promising_child()`` method of the ``MctsNode`` class.
  b. If there is such a child set current node to the child

3.	Return current node

The loop condition in 2. will use the ``is_expanded`` flag of the current node and the ``is_solved`` flag of the state of the current node (see below). 2.a. might not return any child if all the children of the current node were rejected by the tree search (the templates were not applicable). 

Second, the selected leaf node is expanded. This is called the `Expansion` phase in the literature and is used to add new children to a node. The ``expand()`` method of the ``MctsNode`` class takes care of this, but it actually does not instantiate any children nodes. What it does is to use the expansion policy to extract ``RetroReaction`` objects and the probability of each such action. The probability for each action will also be the initial value of each child node. The ``expand()`` method will also set the visitation count for each child to 1. If the ``is_expanded`` flag of the ``MctsNode`` is set or if the ``is_expandable`` flag is not set (see below), the ``expand()`` method will not do anything.

Third, we enter the inner loop of the tree search or the `Rollout` phase, which has the purpose of expanding the tree until we reach a terminal state., i.e. until the ``is_terminal`` flag of the current leaf node is set. The inner loop will execute the following steps:

1.	Retrieve the most promising child of the current leaf node (see above). 
2.	If such a child exists, expand it using the ``expand()`` method and the set current leaf node to this child. 

If 1. does return any child, the ``is_terminal`` flag of the leaf node will have been set and therefore the inner loop will break. Similarly, if the child returned by 1. and set to the current leaf in 2. contains a terminal state, the loop will break. 

Fourth, and finally the algorithm enters the `Backpropagation` phase, which is used to update the value of each node, from the current leaf node all the way to the root. This is done by calling the ``backpropagate()`` method of the ``MctsTreeSearch`` class, which in turn will call the ``backpropagate()`` method of each node on the path between the current leaf and the root.

A few things are worth mentioning about the ``promising_child()`` method of the ``MctsNode`` class. If will select the most promising child by sorting them on the upper confidence bound (UCB) score. The child with the highest score will be selected for instantiation, which means that the ``RetroReaction`` associated with the child will be applied to create new precursors. These precursors will form the state of the new ``MctsNode`` object that is the child. If the application of the reaction failed to produce any precursors, the child value will be set to a large negative value that prevents it from being selected again. The child value will be set to a large negative value also if a filter policy is used in the search and the filter rejects the reaction. Furthermore, ``promising_child()`` will be called recursively until a child can be instantiated (the reaction can be applied). If none of the children can be instantiated the ``is_expanded`` and ``expandable`` flags are updated, and the method returns no child (``None``).

This list explains the different flags of the Node and State objects that are used at various points in the tree search algorithm

================= ============= =========================================================================================== ============================================ ================================================================================= ================================================================================================
Flag              Class         Description                                                                                 Used when                                    Initialized to                                                                    Changed by
================= ============= =========================================================================================== ============================================ ================================================================================= ================================================================================================
``is_expanded``   ``MctsNode``  is True when the node has been expanded                                                     Selection, Expansion                         False when node is created                                                        ``expand()``, sets it to True. ``promising_child()`` sets it to False if no child could be instantiated.
``is_expandable`` ``MctsNode``  is True if the node can be expanded                                                         Rollout (indirectly), Expansion              Set to the opposite of the ``is_terminal`` flag of the state when node is created ``promising_child()`` sets it to False if no child could be instantiated.
``is_terminal()`` ``MctsNode``  is True if either the node is unexpandable or the ``is_terminal`` flag of the state is set  Rollout                                      N/A                                                                               N/A
``is_solved``     ``MctsState`` is True if all precursors in the state is in stock                                          Selection, State init(indirectly)            True or False depending on stock                                                  Never
``is_terminal``   ``MctsState`` is True is ``is_solved`` is True or maximum number of transforms has been reached           Rollout (indirectly), Node init (indirectly) True or False                                                                     Never
================= ============= =========================================================================================== ============================================ ================================================================================= ================================================================================================