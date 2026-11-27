**TomoPy** is an open-source
Python package for tomographic data processing and image reconstruction.

.. image:: https://readthedocs.org/projects/tomopy/badge/?version=latest
   :target: https://readthedocs.org/projects/tomopy/?badge=latest
   :alt: Read the Docs

.. image:: https://dev.azure.com/tomopy/tomopy/_apis/build/status/tomopy.tomopy?branchName=master
   :target: https://dev.azure.com/tomopy/tomopy
   :alt: Azure Pipelines

.. image:: https://coveralls.io/repos/tomopy/tomopy/badge.svg?branch=master
   :target: https://coveralls.io/r/tomopy/tomopy?branch=master
   :alt: Coveralls

.. image:: https://codeclimate.com/github/tomopy/tomopy/badges/gpa.svg
   :target: https://codeclimate.com/github/tomopy/tomopy
   :alt: Code Climate

.. image:: https://anaconda.org/conda-forge/tomopy/badges/downloads.svg
   :target: https://anaconda.org/conda-forge/tomopy
   :alt: Anaconda downloads

.. image:: https://mybinder.org/badge.svg
   :target: https://mybinder.org/v2/gh/tomopy/tomopy/master
   :alt: Use on Binder

|

Features
========

* Image reconstruction algorithms for tomography.
* Various filters, ring removal algorithms, phase retrieval algorithms.
* Forward projection operator for absorption and wave propagation.

Installation
============

This branch provides the ``1.15.3-densq.1`` gridrec fixes as a downloadable
Linux x86-64 package with CUDA 12.9 support and FFTW. With Conda installed,
activate your lftomo environment and run:

.. code-block::

  CONDA_AGGRESSIVE_UPDATE_PACKAGES='' conda install \
    --freeze-installed --override-channels \
    -c https://raw.githubusercontent.com/den-sq/tomopy/1.15.3-densq.1/conda/channel \
    -c conda-forge \
    'libtomo=1.15.3=cuda12.9h1234567_densq1_103' \
    'tomopy=1.15.3=pyhd8ed1ab_2'

Conda downloads the patched native ``libtomo`` from this branch and the
unchanged Python ``tomopy`` wrapper and dependencies from conda-forge. No local
compilation is needed. See the `package instructions <conda/README.rst>`_ for
requirements, the direct download, and the build pin to retain in your environment.


Try Now
=======

Some `Example Jupyter notebooks <https://tomopy.readthedocs.io/en/latest/examples.html>`_ will run in your
browser using `binder <https://mybinder.org/v2/gh/tomopy/tomopy/master?filepath=doc%2Fsource%2Fipynb%2Ftomopy.ipynb>`_.

Contribute
==========

* Bugs Tracker: https://github.com/tomopy/tomopy/projects/1
* Issue Tracker: https://github.com/tomopy/tomopy/issues
* Documentation: https://github.com/tomopy/tomopy/tree/master/doc
* Source Code: https://github.com/tomopy/tomopy/tree/master/tomopy
* Tests: https://github.com/tomopy/tomopy/tree/master/test

License
=======

The project is licensed under the
`BSD-3 <https://github.com/tomopy/tomopy/blob/master/LICENSE.txt>`_ license.
