owls-app
========

Interactive visualization tool for time-series spectra from the Olin Wilson Legacy Survey 
(OWLS).

OWLS is monitoring stellar magnetic activity cycles with optical echelle spectra of 
FGKM stars with the ARC 3.5 m Telescope at Apache Point Observatory (APO).

Getting started
---------------

Installation
^^^^^^^^^^^^

.. code-block:: bash

    python -m pip install owls-app

Run the app
^^^^^^^^^^^

To run the interactive tool, run:

.. code-block:: bash

    owls-app


Multiple sessions
^^^^^^^^^^^^^^^^^

If you already have an instance of `solara` or `owls-app` running, you can
launch a second one by specifying a different port:

.. code-block:: bash

    owls-app --port 1234

Details
-------

When you select a new observation, ``owls-app`` will download and cache the spectrum 
from a collection OWLS FITS files hosted online. The cache is managed using the 
`astropy cache machinery <https://docs.astropy.org/en/stable/utils/data.html>`_, which 
will store the files in a hidden directory in your home directory, ``~/.owls-app/``.
The spectra are stored as ``.fits.gz`` files, and they are 0.63 MB per spectrum. 

To clear the OWLS observations from the cache on your machine, you can call:

.. code-block:: python

    from astropy.utils.data import clear_download_cache
    clear_download_cache(pkgname='owls-app')  # note that the pkgname is important!

Note: if no package name is specified, this command would delete any files cached 
by astropy, and you might not want that.
