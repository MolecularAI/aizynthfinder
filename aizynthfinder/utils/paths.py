""" Module containing routines for returning package paths
"""
import os


def package_path() -> str:
    """Return the path to the package"""
    return os.path.abspath(os.path.dirname(os.path.dirname(__file__)))


def data_path() -> str:
    """Return the path to the ``data`` directory of the package"""
    return os.path.join(package_path(), "data")
