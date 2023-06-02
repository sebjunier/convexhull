#Fichier setup du package convexhull

from distutils.core import setup

setup(name="convexhull",
      version="3.2.0",
      description="Calcul d'enveloppe convexe",
      author="Sebastien Junier",
      # packages created
      packages=['convexhull'],
      #
      install_requires=["numpy", "scipy"])
      # install parameters...
      #license="Apache 2.0")
