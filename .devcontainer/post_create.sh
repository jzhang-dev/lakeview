#!/bin/sh

# Install development dependencies
pip install -v -r dev/requirements.txt

# Install Lakeview for development
pip install -ve .

# Install jupyterlab-horizon-theme
pip install git+https://github.com/jzhang-dev/jupyterlab-horizon-theme

# Build docs
sphinx-build -W --keep-going -E -b html docs/ docs/_build/html
