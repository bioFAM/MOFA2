import sys
from setuptools import setup
from setuptools import find_packages

def setup_package():
  install_requires = ['pandas', 'scipy', 'numpy', 'sklearn', 'argparse', 'h5py']
  metadata = dict(
      name = 'mofapy2',
      version = '0.1',
      description = 'Multi-Omics Factor Analysis v2, a statistical framework for the integration of multi-group and multi-omics data',
      url = 'http://github.com/PMBio/biofam',
      author = 'Ricard Argelaguet',
      author_email = 'ricard.argelaguet@gmail.com',
      license = 'LGPL-3.0',
      packages = find_packages(),
      install_requires = install_requires
    )

  setup(**metadata)


if __name__ == '__main__':
  if sys.version_info < (2,7):
    sys.exit('Sorry, Python < 2.7 is not supported')
    
  setup_package()