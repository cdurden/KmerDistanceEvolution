from setuptools import setup, find_packages

requires = [
#    'Bio',
#    'scipy',
    'sqlalchemy',
    'numpy',
    'biopython',
]

setup(name='Coalescent',
      version='1.0',
      install_requires=requires,
      package_dir = {'Coalescent': '.'},
      packages=find_packages('..'),
)

