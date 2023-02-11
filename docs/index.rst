.. Lakeview documentation master file, created by
   sphinx-quickstart on Fri Dec 30 17:31:55 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. figure:: gallery/output/logo.svg
   :align: left
   :width: 600
   :alt: Lakeview logo

Welcome to Lakeview's documentation!
====================================


Lakeview is a Python 3 library for creating publication-quality `IGV <https://software.broadinstitute.org/software/igv/>`_-style genomic visualizations. Lakeview is based on `Matplotlib <https://matplotlib.org/>`_. 



Features
------------

- **Improved clarity**. Lakeview inherits the familiar and intuitive visual style of `IGV <https://software.broadinstitute.org/software/igv/>`_, with a clear layout designed for publication and presentation. 
- **Programmable plotting**. Multiple files and genomic regions can be visualized automatically through a Pythonic interface inspired by `Seaborn <https://seaborn.pydata.org/>`_ and `Pandas <https://pandas.pydata.org/>`_.
- **Support for remote data**. Genomic data are often stored in remote servers without display devices. With Lakeview, you can plot remotely and view the output figures locally. Lakeview works well with `JupyterLab <https://jupyterlab.readthedocs.io/en/stable/>`_ to streamline this workflow. 
- **Transparency and reproduciblity**. Figures are plotted transparently and annotated explicitly. The input data and the plotting code contain all the information needed to reproduce the figure. 
- **Customizable layouts**. Lakeview supports many layouts implemented in `IGV <https://software.broadinstitute.org/software/igv/>`_, while allowing the user to define custom rules for ordering, groupping, and coloring each segment. Advanced customization is possible via the `Matplotlib <https://matplotlib.org/>`_ API.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tutorials/index
   gallery/index
   api/index


