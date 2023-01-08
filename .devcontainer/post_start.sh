#!/bin/sh

# Override JupyterLab settings
mkdir -p /opt/conda/share/jupyter/lab/settings
cp dev/jupyterlab_settings/overrides.json /opt/conda/share/jupyter/lab/settings/overrides.json

# Start JupyterLab server
jupyter lab --no-browser --ip=0.0.0.0 --port=9999 --NotebookApp.token='' --NotebookApp.password=''
