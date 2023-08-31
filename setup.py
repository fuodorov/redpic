from setuptools import find_packages, setup

from redpic.core import config

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements/base.txt", "r", encoding="utf-8") as fh:
    required = fh.read().splitlines()

setup(
    name=config.NAME,
    version=config.VERSION,
    author=config.AUTHOR,
    author_email=config.AUTHOR_EMAIL,
    description=config.DOC,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=config.URL,
    packages=find_packages(),
    install_requires=required,
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    project_urls={
        "Bug Reports": "https://github.com/fuodorov/redpic/issues",
        "PyPi": "https://pypi.org/project/redpic/",
        "Documentation": "https://fuodorov.github.io/redpic/",
        "Source": "https://github.com/fuodorov/redpic/",
    },
)
