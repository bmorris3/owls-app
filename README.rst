owls-app
========

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

.. image:: http://img.shields.io/badge/powered%20by-jdaviz-336699.svg?style=flat
    :target: https://github.com/spacetelescope/jdaviz/
    :alt: Powered by jdaviz

.. image:: http://img.shields.io/badge/powered%20by-lightkurve-336699.svg?style=flat
    :target: https://github.com/spacetelescope/lcviz/
    :alt: Powered by lcviz

.. image:: http://img.shields.io/badge/powered%20by-specutils-ff9722.svg?style=flat
    :target: https://github.com/astropy/specutils/
    :alt: Powered by specutils


Interactive visualization tool for time-series spectra from the Olin Wilson Legacy Survey 
(OWLS).

OWLS is monitoring stellar magnetic activity cycles with optical echelle spectra of 
FGKM stars with the ARC 3.5 m Telescope at Apache Point Observatory (APO).


.. image:: https://github.com/bmorris3/owls-app/blob/main/docs/owls_demo.png?raw=true
    :alt: screenshot of owls-app


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

What is this spectrum viewer?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Meet **Specviz** (`documentation <https://jdaviz.readthedocs.io/en/stable/specviz/index.html>`__,
`source <https://github.com/spacetelescope/jdaviz/>`__) 
part of the ``jdaviz`` interactive data visualization and analysis package 
developed at the Space Telescope Science Institute. Specviz supports quite
advanced workflows that we won't summarize here -- check out their docs at
the link above.

What is this time series viewer?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Meet **LCviz** (`documentation <https://lcviz.readthedocs.io/>`__, 
`source <https://github.com/spacetelescope/lcviz>`__),  a light curve visualization
and analysis tool built on ``jdaviz``, developed at the Space Telescope Science 
Institute. LCviz supports quite advanced workflows that we won't summarize
here -- check out their docs at the link above.


How are the echelle orders normalized?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The app will attempt to remove the continuum from each echelle order by 
dividing out a fifth order polynomial fit to the corresponding order of 
an ARCES spectrum of the hot subdwarf standard
`HZ 44 <https://simbad.cds.unistra.fr/simbad/sim-id?Ident=HZ+44>`__ to roughly
remove the blaze function, and then dividing by the order's maximum flux.
This normalization method is most imprecise near strong lines in hot star 
atmospheres. The flux units of the spectrum viewer are thus relative flux.


Can I download these spectra?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yes! We're submitting all spectra to MAST as HLSPs, though it may take some time
before they're available on MAST. `Reach out to Brett <mailto:morrisbrettm@gmail.com>`__
if you'd like them sooner.


Are the spectra saved locally?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When you select a new observation, ``owls-app`` will download and cache the spectrum 
from a collection OWLS FITS files hosted online. The cache is managed using the 
`astropy cache machinery <https://docs.astropy.org/en/stable/utils/data.html>`__, which 
will store the files in a hidden directory in your home directory, ``~/.owls-app/``.
The spectra are stored as ``.fits.gz`` files, and they are 0.63 MB per spectrum. 

To clear the OWLS observations from the cache on your machine, you can call:

.. code-block:: python

    from astropy.utils.data import clear_download_cache
    clear_download_cache(pkgname='owls-app')  # note that the pkgname is important!

Note: if no package name is specified, this command would delete any files cached 
by astropy, and you might not want that.


