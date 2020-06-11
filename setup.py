from setuptools import setup, find_packages

with open("requirements.txt", "r") as fileobj:
    requirements = [req.strip() for req in fileobj.readlines()]

setup(
    name="aizynthfinder",
    version="1.0.0",
    description="Retrosynthetic route finding using neural network guided Monte-Carlo tree search.",
    author="Molecular AI group",
    author_email="samuel.genheden@astrazeneca.com",
    license="MIT",
    packages=find_packages(exclude=("tests",)),
    install_requires=requirements,
    package_data={"aizynthfinder": ["data/*.yml"]},
    entry_points={
        "console_scripts": [
            "aizynthapp = aizynthfinder.interfaces.aizynthapp:main",
            "aizynthcli = aizynthfinder.interfaces.aizynthcli:main",
            "smiles2stock = aizynthfinder.tools.make_stock:main",
            "preprocess_rollout = aizynthfinder.training.preprocess_rollout:main",
            "aizynth_training = aizynthfinder.training.training:main",
        ]
    },
)
