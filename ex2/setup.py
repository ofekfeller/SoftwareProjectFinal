from setuptools import setup, Extension

"""
A minimalist setup is shown.
"""


setup(name='mykmeanssp',
      version='1.0',
      description='moth**fuc**ingclass',
      ext_modules=[Extension('mykmeanssp', sources=['kmeans.c'])])
