""" Module containing routine to dynamically load a class from a specification """
from __future__ import annotations
import importlib
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any


def load_dynamic_class(
    name_spec: str, default_module: str = None, exception_cls: Any = ValueError
) -> Any:
    """
    Load an object from a dynamic specification.

    The specification can be either:
        ClassName, in-case the module name is taken from the `default_module` argument
    or
        package_name.module_name.ClassName, in-case the module is taken as `package_name.module_name`

    :param name_spec: the class specification
    :param default_module: the default module
    :param exception_cls: the exception class to raise on exception
    :return: the loaded class
    """
    if "." not in name_spec:
        name = name_spec
        if not default_module:
            raise exception_cls(
                "Must provide default_module argument if not given in name_spec"
            )
        module_name = default_module
    else:
        module_name, name = name_spec.rsplit(".", maxsplit=1)

    try:
        loaded_module = importlib.import_module(module_name)
    except ImportError:
        raise exception_cls(f"Unable to load module: {module_name}")

    if not hasattr(loaded_module, name):
        raise exception_cls(
            f"Module ({module_name}) does not have a class called {name}"
        )

    return getattr(loaded_module, name)
