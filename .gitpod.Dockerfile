FROM gitpod/workspace-full:2022-09-11-15-11-40

RUN pip install Bio black ipywidgets jupyterlab matplotlib ncbi-datasets-pylib pysam pytest sphinx==4.5.0 sphinx-book-theme==0.3.3
RUN brew install samtools