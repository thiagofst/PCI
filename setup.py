from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('license.txt') as f:
    license = f.read()

setup(
    name='pci',
    version='1.0.0',
    description='Propagation of Confidence Intervals (PCI). Python module for dealing and propagate measurments and its assymetric confidence levels',
    long_description=readme,
    author='Muryel Guolo Pereira',
    author_email='muryel@astro.ufsc.br',
    url='https://github.com/muryelgp/PCI',
    download_url = 'https://github.com/muryelgp/PCI/archive/master.zip',
    license=license,
    keywords = ['statistics', 'uncertanties', 'confidence-level'],
    packages=['pci'],
    install_requires=[
          'numpy',
          'matplotlib',
          'scipy',
      ],
)