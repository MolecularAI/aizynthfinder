""" Module containing custom exception classes
"""


class CostException(Exception):
    """Exception raised by molecule cost classes"""


class ExternalModelAPIError(Exception):
    """ Custom error type to signal failure in External model"""


class MoleculeException(Exception):
    """An exception that is raised by molecule class"""


class NodeUnexpectedBehaviourException(Exception):
    """Exception that is raised if the tree search is behaving unexpectedly."""


class PolicyException(Exception):
    """An exception raised by policy classes"""


class RejectionException(Exception):
    """ An exception raised if a retro action should be rejected"""


class ScorerException(Exception):
    """Exception raised by scoring classes"""


class StockException(Exception):
    """An exception raised by stock classes"""


class TreeAnalysisException(Exception):
    """Exception raised when analysing trees"""
