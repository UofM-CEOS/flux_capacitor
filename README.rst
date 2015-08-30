Routines and tools to process flux (eddy covariance) data as collected
by our group.

Basic usage
===========

The package currently has 2 sub-packages:

-  ``eddycov``
-  ``underway``

=eddycov= is the eddy covariance package, and ``underway`` is for the
calculation of pCO2 from the underway system.

Both packages can be imported at once by:

.. code:: python

      import fluxer

whereby each sub-package's namespace is accessible via
``fluxer.eddycov.*`` or ``fluxer.underway.*``.

However, each package can be be imported individually, if needed:

.. code:: python

      import fluxer.eddycov as eddycov

Or

.. code:: python

      import fluxer.underway as underway

thus avoiding having to use the ``fluxer.`` prefix.

Project configuration
=====================

The easiest way to use the packages is to set up a configuration for any
given project. The source for ``fluxer`` includes the default
configuration settings required by each sub-package (under the
``config/`` directory). These settings are specified in a ``*.cfg`` file
(syntax instructions are given in the default files).

=eddycov= package
=================

The main interface for this package is the two functions: ``main`` and
``flux_period``:

.. code:: python

      # The main() function takes a configuration file and runs the analyses
      from fluxer import eddycov
      eddycov.main("ec_config.cfg")

.. code:: python

      # The flux_period() function takes an input data file and a *parsed*
      # configuration file and runs the analysis for it
      from fluxer import eddycov
      from fluxer.flux_config import parse_config
      config = parse_config("ec_config.cfg")
      eddycov.flux_period("YYYYMMDD_100000_10hz.csv", config)

A third alternative is offered to allow execution from a shell command
line:

.. code:: shell

      db_flux.py ec_config.cfg

However, this requires ensuring that the shebang (first line of the
script/module) on ``db_flux.py`` is appropriate for the system running
it.

=underway= package
==================

Consistent with ``eddycov``, this package offers two functions: ``main``
and ``underway_pCO2``.

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

