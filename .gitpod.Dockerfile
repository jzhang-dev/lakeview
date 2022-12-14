FROM gitpod/workspace-full:2022-09-11-15-11-40

# RUN brew install samtools
RUN pip install Bio black ipywidgets jupyterlab matplotlib ncbi-datasets-pylib pysam pytest
RUN pip install git+https://github.com/jzhang-dev/jupyterlab-horizon-theme