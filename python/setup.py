import configparser
import os.path

import pybind11
import setuptools


def get_metadata(key):
    with open(os.path.join("..", "metadata", key+".txt"), "r", encoding="utf8") as f:
        return f.read().strip()


if __name__ == "__main__":
    metadata = {
        "name": "hd98",
        "version": get_metadata("version"),
        "author": get_metadata("author"),
        "author_email": "email",
        "description": get_metadata("description"),
        "url": get_metadata("repository"),
    }

    with open(os.path.join("..", "README.md"), "r") as f:
        metadata["long_description"] = f.read()

    config = configparser.ConfigParser()
    config.read("setup.cfg")
    hd98_include_dir = config["hd98"].get("include_dir", "")
    hd98_library_dir = config["hd98"].get("library_dir", "")

    hd98 = setuptools.Extension(
        "hd98.hd98",
        include_dirs=[pybind11.get_include(),
                      hd98_include_dir],
        sources=[os.path.join("hd98",
                              "hd98.cpp")],
        libraries=["hd98"],
        library_dirs=[hd98_library_dir],
    )

    setuptools.setup(
        long_description_content_type="text/markdown",
        packages=setuptools.find_packages(),
        ext_modules=[hd98],
        **metadata
    )
