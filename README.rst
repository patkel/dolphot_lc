==========
Dolphot-LC
==========

|PyPI Status| |RTD Status|

Dolphot-LC is an automated Hubble Space Telescope (HST) data pipeline based on the popular Dolphot analysis package. This package allows those interested in astronomy to create lightcurves and difference images from HST fits images. 

Installation
-------

Install dolphot_lc with pip.


``pip install dolphot_lc``

Detailed installation instructions can be found on our ReadTheDocs.    

Getting Started
-------

After following the full installation, test success by running the following. 

``python3 dolphot_lc_test.py``


Documentation
-------

Installation and development instructions can be found on our Read the Docs and Jupyter Notebook. The Jupyter Notebook allows the user to run their own fits images through the pipeline and generate results. Dolphot-LC requires a coadded template image and science images that are already aligned to template; our testing procedures have example images already in the github. 

Read the Docs ---- https://dolphot-lc.readthedocs.io/en/latest/

Jupyter Notebook ---- https://nbviewer.org/gist/whit5224/287af111f44bf83a23eaaf19a5121c75

Contributing
-------

Contributions are always welcome! See `contributing.md` for ways to get started.

Please adhere to this project's `code of conduct`.

License
-------

Dolphot-LC is licensed under a 3-clause BSD style license - see the
`LICENSE.rst <LICENSE.rst>`_ file.

.. |PyPI Status| image:: https://img.shields.io/pypi/v/dolphot_lc.svg
    :target: https://pypi.org/project/dolphot-lc/
    :alt: Dolphot-LC's PyPI Status

.. |RTD Status| image:: https://readthedocs.org/projects/dolphot-lc/badge/?version=latest
    :target: https://dolphot-lc.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
