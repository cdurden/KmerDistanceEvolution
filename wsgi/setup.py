from setuptools import setup, find_packages

requires = [
    'sqlalchemy',
    'venusian',
    'pyramid_chameleon',
    'pyramid',
    'pyramid_tm',
    'waitress',
    'zope.sqlalchemy',
    'matplotlib',
    'Coalescent'
]

setup(name='kmers',
      install_requires=requires,
      packages=find_packages(),
      package_data={
          '': ['templates/*.pt']
          },
      entry_points="""\
      [paste.app_factory]
      main = kmers:main
      [console_scripts]
      initdb = kmers.initdb:main
      """,
)
