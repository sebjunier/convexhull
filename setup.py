#Fichier setup du package convexhull

from setuptools import setup, find_packages

setup(
name="convexhull",
      version="4.0.1",
      description="Convexhull calculation",
      author="Sebastien Junier",
      # packages created
      packages=find_packages(),
      # dependances
      install_requires=["numpy"],
      # install parameters...
      #license="Apache 2.0",
      )
