import pytest

from aizynthfinder.context.collection import ContextCollection


class StringCollection(ContextCollection):
    def load(self, key, value):
        self._items[key] = value

    def load_from_config(self, config):
        for key, value in config.items():
            self.load(key, value)


class SingleStringCollection(StringCollection):
    _single_selection = True


def test_empty_collection():
    collection = StringCollection()

    assert len(collection) == 0
    assert collection.selection == []


def test_add_single_item():
    collection = StringCollection()

    collection.load("key1", "value1")

    assert len(collection) == 1
    assert collection.items == ["key1"]
    assert collection.selection == []


def test_get_item():
    collection = StringCollection()
    collection.load("key1", "value1")

    assert collection["key1"] == "value1"

    with pytest.raises(KeyError):
        collection["key2"]


def test_del_item():
    collection = StringCollection()
    collection.load("key1", "value1")

    del collection["key1"]

    assert len(collection) == 0

    with pytest.raises(KeyError):
        del collection["key2"]


def test_select_single_item():
    collection = StringCollection()
    collection.load("key1", "value1")

    collection.selection = "key1"

    assert collection.selection == ["key1"]

    with pytest.raises(KeyError):
        collection.selection = "key2"

    collection.load("key2", "value2")

    collection.selection = "key2"

    assert collection.selection == ["key2"]


def test_select_append():
    collection = StringCollection()
    collection.load("key1", "value1")
    collection.load("key2", "value2")

    collection.selection = "key1"

    assert collection.selection == ["key1"]

    collection.select("key2", append=True)

    assert collection.selection == ["key1", "key2"]


def test_select_all():
    collection = StringCollection()

    collection.select_all()

    assert collection.selection == []

    collection.load("key1", "value1")
    collection.load("key2", "value2")

    collection.select_all()

    assert collection.selection == ["key1", "key2"]


def test_select_first():
    collection = StringCollection()

    collection.select_first()

    assert collection.selection == []

    collection.load("key1", "value1")
    collection.load("key2", "value2")

    collection.select_first()

    assert collection.selection == ["key1"]


def test_select_overwrite():
    collection = StringCollection()
    collection.load("key1", "value1")
    collection.load("key2", "value2")

    collection.selection = "key1"

    assert collection.selection == ["key1"]

    collection.selection = ["key2"]

    assert collection.selection == ["key2"]


def test_deselect_all():
    collection = StringCollection()
    collection.load("key1", "value1")
    collection.load("key2", "value2")
    collection.selection = ("key1", "key2")

    assert len(collection) == 2
    assert collection.selection == ["key1", "key2"]

    collection.deselect()

    assert len(collection.selection) == 0


def test_deselect_one():
    collection = StringCollection()
    collection.load("key1", "value1")
    collection.load("key2", "value2")
    collection.selection = ("key1", "key2")

    collection.deselect("key1")

    assert collection.selection == ["key2"]

    with pytest.raises(KeyError):
        collection.deselect("key1")


def test_empty_single_collection():
    collection = SingleStringCollection()

    assert len(collection) == 0
    assert collection.selection is None


def test_select_single_collection():
    collection = SingleStringCollection()
    collection.load("key1", "value1")
    collection.load("key2", "value2")

    collection.selection = "key1"

    assert collection.selection == "key1"

    collection.selection = "key2"

    assert collection.selection == "key2"


def test_select_multiple_single_collection():
    collection = SingleStringCollection()
    collection.load("key1", "value1")
    collection.load("key2", "value2")

    collection.selection = ["key1"]

    assert collection.selection == "key1"

    with pytest.raises(ValueError):
        collection.selection = ["key1", "key2"]
