from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()

# FLUXER_CONFIGS = ["fluxer/config/flux_default.cfg",
#                   "fluxer/config/underway_default.cfg"]
setup(name="fluxer",
      version="0.1.0",
      description=("Tools to process flux (eddy covariance) data"
                   "collected by CEOS."),
      long_description=readme(),
      author="Sebastian Luque",
      author_email="sebastian.luque@umanitoba.ca",
      url="https://github.com/UofM-CEOS/flux_capacitor",
      packages=["fluxer", "fluxer.eddycov", "fluxer.underway"],
      scripts=["pCO2/pCO2_rosette.py", "pCO2/subset_bottles.awk",
               "pCO2/pCO2_bottle_match.awk"],
      test_suite="fluxer.tests",
      entry_points={
          "console_scripts": ["get_fluxes = fluxer.eddycov.__main__:main"]
      }
      # package_data={'fluxer': FLUXER_CONFIGS},
      # data_files=[('config', ['config/flux_default.cfg',
      #                         'config/underway_default.cfg'])],
      )
