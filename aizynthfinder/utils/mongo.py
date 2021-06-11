""" Module containing routines to obtain a MongoClient instance
"""
from urllib.parse import urlencode

from pymongo import MongoClient

from aizynthfinder.utils.logging import logger

_CLIENT = None


def get_mongo_client(
    host: str = "localhost",
    port: int = 27017,
    user: str = None,
    password: str = None,
    tls_certs_path: str = "",
) -> MongoClient:
    """
    A helper function to create and reuse MongoClient

    The client is only setup once. Therefore if this function is called a second
    time with different parameters, it would still return the first client.

    :param host: the host
    :param port: the host port
    :param user: username, defaults to None
    :param password: password, defaults to None
    :param tls_certs_path: the path to TLS certificates if to be used, defaults to ""
    :raises ValueError: if host and port is not given first time
    :return: the MongoDB client
    """
    global _CLIENT
    if _CLIENT is None:
        params = {}
        if tls_certs_path:
            params.update({"ssl": "true", "ssl_ca_certs": tls_certs_path})
        cred_str = f"{user}:{password}@" if password else ""
        uri = f"mongodb://{cred_str}{host}:{port}/?{urlencode(params)}"
        logger().debug(f"Connecting to MongoDB on {host}:{port}")
        _CLIENT = MongoClient(uri)
    return _CLIENT
