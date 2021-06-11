""" Module with script to download public data
"""
import argparse
import os
import sys

import requests
import tqdm


FILES_TO_DOWNLOAD = {
    "policy_model": {
        "filename": "uspto_model.hdf5",
        "url": "https://ndownloader.figshare.com/files/23086454",
    },
    "template_file": {
        "filename": "uspto_templates.hdf5",
        "url": "https://ndownloader.figshare.com/files/23086457",
    },
    "stock": {
        "filename": "zinc_stock.hdf5",
        "url": "https://ndownloader.figshare.com/files/23086469",
    },
    "filter_policy": {
        "filename": "uspto_filter_model.hdf5",
        "url": "https://ndownloader.figshare.com/files/25584743",
    },
}

YAML_TEMPLATE = """policy:
  files:
    uspto:
      - {}
      - {}
filter:
  files:
    uspto: {}
stock:
  files:
    zinc: {}
"""


def _download_file(url: str, filename: str) -> None:
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        total_size = int(response.headers.get("content-length", 0))
        pbar = tqdm.tqdm(
            total=total_size, desc=os.path.basename(filename), unit="B", unit_scale=True
        )
        with open(filename, "wb") as fileobj:
            for chunk in response.iter_content(chunk_size=1024):
                fileobj.write(chunk)
                pbar.update(len(chunk))
        pbar.close()


def main() -> None:
    """Entry-point for CLI"""
    parser = argparse.ArgumentParser("download_public_data")
    parser.add_argument(
        "path",
        default=".",
        help="the path download the files",
    )
    path = parser.parse_args().path

    try:
        for filespec in FILES_TO_DOWNLOAD.values():
            _download_file(filespec["url"], os.path.join(path, filespec["filename"]))
    except requests.HTTPError as err:
        print(f"Download failed with message {str(err)}")
        sys.exit(1)

    with open(os.path.join(path, "config.yml"), "w") as fileobj:
        path = os.path.abspath(path)
        fileobj.write(
            YAML_TEMPLATE.format(
                os.path.join(path, FILES_TO_DOWNLOAD["policy_model"]["filename"]),
                os.path.join(path, FILES_TO_DOWNLOAD["template_file"]["filename"]),
                os.path.join(path, FILES_TO_DOWNLOAD["filter_policy"]["filename"]),
                os.path.join(path, FILES_TO_DOWNLOAD["stock"]["filename"]),
            )
        )
    print("Configuration file written to config.yml")


if __name__ == "__main__":
    main()
