from setuptools import setup
from setuptools import find_packages

def setup_package():
  install_requires = []
  metadata = dict(
      name = 'new_model_name',
      version = '0.1',
      description = 'Factor Analysis model for ',
      #long_description=read('README.rst'),
      packages = find_packages(),
    )

  setup(**metadata)

if __name__ == '__main__':
  setup_package()
