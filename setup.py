from setuptools import setup
import os

defs = os.path.join('src')

long_description = open('./README.md', 'r').read()

setup(name='MCSPS',
      version='0.1',
      description='A Python tool for performing site permutation search on fixed atomic geometries.',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Frank T. Cerasoli',
      author_email='ftcerasoli@ucdavis.edu',
      platforms='Unix',
      url='https://github.com/nanotheorygroup/MC-SPS',
      packages=['MCSPS'],
      package_dir={'MCSPS':'src'},
      install_requires=['ase', 'numpy'],
      extra_requires=['matplotlib', 'megnet'],
      python_requires='>=3.7'
)
