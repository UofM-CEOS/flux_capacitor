from distutils.core import setup

setup(name='fluxer',
      version='0.1.0',
      description='Tools to process flux (eddy covariance) data collected by CEOS.',
      author='Sebastian Luque',
      author_email='sebastian.luque@umanitoba.ca',
      url='https://github.com/UofM-CEOS/flux_capacitor',
      packages=['fluxer', 'fluxer.eddycov', 'fluxer.underway'],
      # package_data={'flux': ['flux/config/flux_default.cfg']},
      # data_files=[('config', ['config/flux_default.cfg',
      #                         'config/underway_default.cfg'])],
      )
