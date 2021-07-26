import pytest

from aizynthfinder.utils.loading import load_dynamic_class


class MyException(Exception):
    pass


def test_simple():
    cls = load_dynamic_class("aizynthfinder.reactiontree.ReactionTree")

    assert cls.__name__ == "ReactionTree"


def test_default_module():
    cls = load_dynamic_class(
        "ReactionTree", default_module="aizynthfinder.reactiontree"
    )

    assert cls.__name__ == "ReactionTree"


def test_no_default_module():
    with pytest.raises(ValueError, match="default_module"):
        load_dynamic_class("ReactionTree")

    with pytest.raises(MyException, match="default_module"):
        load_dynamic_class("ReactionTree", exception_cls=MyException)


def test_incorrect_module():
    bad_module = "aizynthfinder.rt."
    with pytest.raises(ValueError, match=bad_module):
        load_dynamic_class(f"{bad_module}.ReactionTree")

    with pytest.raises(MyException, match=bad_module):
        load_dynamic_class(f"{bad_module}.ReactionTree", exception_cls=MyException)


def test_incorrect_class():
    bad_class = "ReactionTreee"
    with pytest.raises(ValueError, match=bad_class):
        load_dynamic_class(f"aizynthfinder.reactiontree.{bad_class}")

    with pytest.raises(MyException, match=bad_class):
        load_dynamic_class(
            f"aizynthfinder.reactiontree.{bad_class}", exception_cls=MyException
        )
