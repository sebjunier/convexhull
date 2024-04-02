#Fichier setup du package convexhull

from setuptools import setup, find_packages

setup(
name="convexhull",
      version="3.2.1",
      description="Calcul d'enveloppe convexe",
      author="Sebastien Junier",
      # packages created
      packages=find_packages(),
      # dependances
      install_requires=["numpy", "scipy"],
      # install parameters...
      #license="Apache 2.0",
      )
