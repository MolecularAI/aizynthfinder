""" Module containing a class that is the base class for all collection classes (stock, policies, scorers)
"""
from __future__ import annotations
import abc
from typing import TYPE_CHECKING

from aizynthfinder.utils.logging import logger

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import StrDict
    from aizynthfinder.utils.type_utils import Any, List, Union


class ContextCollection(abc.ABC):
    """
    Abstract base class for a collection of items
    that can be loaded and then (de-)selected.

     One can obtain individual items with:

    .. code-block::

        an_item = collection["key"]

    And delete items with

    .. code-block::

        del collection["key"]


    """

    _single_selection = False
    _collection_name = "collection"

    def __init__(self):
        self._items: StrDict = {}
        self._selection: List[str] = []
        self._logger = logger()

    def __delitem__(self, key: str) -> None:
        if key not in self._items:
            raise KeyError(
                f"{self._collection_name.capitalize()} with name {key} not loaded."
            )
        del self._items[key]

    def __getitem__(self, key: str) -> Any:
        if key not in self._items:
            raise KeyError(
                f"{self._collection_name.capitalize()} with name {key} not loaded."
            )
        return self._items[key]

    def __len__(self) -> int:
        return len(self._items)

    @property
    def items(self) -> List[str]:
        """The available item keys"""
        return list(self._items.keys())

    @property
    def selection(self) -> Union[List[str], str, None]:
        """The keys of the selected item(s)"""
        if self._single_selection:
            return self._selection[0] if self._selection else None
        return self._selection

    @selection.setter
    def selection(self, value: str) -> None:
        self.select(value)

    def deselect(self, key: str = None) -> None:
        """
        Deselect one or all items

        If no key is passed, all items will be deselected.

        :param key: the key of the item to deselect, defaults to None
        :raises KeyError: if the key is not among the selected ones
        """
        if not key:
            self._selection = []
            return

        if key not in self._selection:
            raise KeyError(f"Cannot deselect {key} because it is not selected")
        self._selection.remove(key)

    @abc.abstractmethod
    def load(self, *_: Any) -> None:
        """Load an item. Needs to be implemented by a sub-class"""

    @abc.abstractmethod
    def load_from_config(self, **config: Any) -> None:
        """Load items from a configuration. Needs to be implemented by a sub-class"""

    def select(self, value: Union[str, List[str]], append: bool = False) -> None:
        """
        Select one or more items.

        If this is a single selection collection, only a single value is accepted.
        If this is a multiple selection collection it will overwrite the selection completely,
        unless ``append`` is True and a single key is given.

        :param value: the key or keys of the item(s) to select
        :param append: if True will append single keys to existing selection
        :raises ValueError: if this a single collection and value is multiple keys
        :raises KeyError: if at least one of the keys are not corresponding to a loaded item
        """
        if self._single_selection and not isinstance(value, str) and len(value) > 1:
            raise ValueError(f"Cannot select more than one {self._collection_name}")

        keys = [value] if isinstance(value, str) else value

        for key in keys:
            if key not in self._items:
                raise KeyError(
                    f"Invalid key specified {key} when selecting {self._collection_name}"
                )

        if self._single_selection:
            self._selection = [keys[0]]
        elif isinstance(value, str) and append:
            self._selection.append(value)
        else:
            self._selection = list(keys)

        self._logger.info(f"Selected as {self._collection_name}: {', '.join(keys)}")

    def select_all(self) -> None:
        """Select all loaded items"""
        if self.items:
            self.select(self.items)

    def select_first(self) -> None:
        """Select the first loaded item"""
        if self.items:
            self.select(self.items[0])

    def select_last(self) -> None:
        """Select the last loaded item"""
        if self.items:
            self.select(self.items[-1])
