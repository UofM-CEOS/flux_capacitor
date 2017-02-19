################
 fluxer package
################


Routines and tools to process flux (eddy covariance) data as collected by
our group.


Basic usage
===========

The package currently has 2 sub-packages:

-  ``eddycov``
-  ``underway``

``eddycov`` is the eddy covariance package, and ``underway`` is a package
for calculation of pCO2 from an underway system.  Each package can be
imported individually as usual, e.g.:

.. code:: python

    import fluxer.eddycov as eddycov

Or

.. code:: python

    import fluxer.underway as underway


Project configuration
=====================

The easiest way to use the packages is to set up a configuration for any
given project. The source for ``fluxer`` includes the default configuration
settings required by each sub-package (under the ``config/``
directory). These settings are specified in a ``*.cfg`` file (syntax
instructions are given in the default files).


``eddycov`` package
===================

The main interface for this package is the function ``main``:

.. code:: python

    # The main() function takes a configuration file and runs the analyses
    import fluxer.eddycov as ec
    ec.main("ec_config.cfg")

and there is also a command-line utility:

.. code:: shell

    get_fluxes ec_config.cfg


``underway`` package
====================

This package offers two functions: ``main`` and ``underway_pCO2``.

.. code:: python

    # The main() function takes a configuration file and runs the analyses
    from fluxer import underway
    underway.main("ec_config.cfg")

.. code:: python

    # The underway_pCO2() function takes an input data file and a *parsed*
    # configuration file and runs the analysis for it
    from fluxer import underway
    from fluxer.flux_config import parse_config
    config = parse_config("uw_config.cfg")
    eddycov.underway_pCO2("YYYYMMDD_100000_20min.csv", config)

It is also possible to perform the analysis from the shell command line:

.. code:: shell

    underway.py ec_config.cfg


###################
 API documentation
###################

.. toctree::
   :maxdepth: 4

   fluxer


####################
 Indices and tables
####################


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
