FROM gitpod/workspace-full:2022-09-11-15-11-40

# RUN brew install samtools
RUN pip install Bio black jupyterlab ncbi-datasets-pylib pytest
RUN pip install git+https://github.com/jzhang-dev/jupyterlab-horizon-theme