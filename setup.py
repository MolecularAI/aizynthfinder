from setuptools import setup, find_packages

with open("requirements.txt", "r") as fileobj:
    requirements = [req.strip() for req in fileobj.readlines()]

setup(
    name="aizynthfinder",
    version="2.2.1",
    description="Retrosynthetic route finding using neural network guided Monte-Carlo tree search.",
    author="Molecular AI group",
    author_email="samuel.genheden@astrazeneca.com",
    license="MIT",
    packages=find_packages(exclude=("tests",)),
    install_requires=requirements,
    package_data={"aizynthfinder": ["data/*.yml", "data/templates/*"]},
    entry_points={
        "console_scripts": [
            "aizynthapp = aizynthfinder.interfaces.aizynthapp:main",
            "aizynthcli = aizynthfinder.interfaces.aizynthcli:main",
            "smiles2stock = aizynthfinder.tools.make_stock:main",
            "cat_aizynth_output = aizynthfinder.tools.cat_output:main",
            "preprocess_rollout = aizynthfinder.training.preprocess_expansion:main",
            "preprocess_expansion = aizynthfinder.training.preprocess_expansion:main",
            "preprocess_filter = aizynthfinder.training.preprocess_filter:main",
            "preprocess_recommender = aizynthfinder.training.preprocess_recommender:main",
            "make_false_products = aizynthfinder.training.make_false_products:main",
            "aizynth_training = aizynthfinder.training.training:main",
            "download_public_data = aizynthfinder.tools.download_public_data:main",
        ]
    },
)
