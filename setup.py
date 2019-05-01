from setuptools import setup, find_packages

requires = [
#    'Bio',
    'sqlalchemy',
    'numpy',
    'scipy',
    'biopython',
]

setup(name='KmerDistanceEvolution',
      version='1.0',
      install_requires=requires,
#      package_dir = {'KmerDistanceEvolution': '.'},
      packages=find_packages('.'),
)

