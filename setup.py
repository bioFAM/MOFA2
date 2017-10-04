from setuptools import setup
from setuptools import find_packages

def setup_package():
  install_requires = ['pandas', 'scipy', 'numpy', 'sklearn', 'argparse', 'h5py']
  # console_scripts = [ 'biofam=biofam.build_model.init_model:entry_point'],
  metadata = dict(
      name = 'biofam',
      version = '0.1',
      description = 'Bio Factor Analysis Models',
      #long_description=read('README.rst'),
      url = 'http://github.com/PMBio/biofam',
      author = 'Ricard Argelaguet, Damien Arnol',
      author_email = 'ricard.argelaguet@gmail.com',
      license = 'MIT',
      packages = find_packages(),
      install_requires = install_requires,
      #entry_points = {'console_scripts': console_scripts}
    )

  setup(**metadata)

if __name__ == '__main__':
  setup_package()
