# convexhull

The `convexhull` package computes the convex hull in N dimensions from the energies of a set of phases. It is based on the Quickhull algorithm [1].

The package can be used to:

* identify the phases that lie on the convex hull;
* calculate the energy above the convex hull for each phase;
* determine the phase decomposition at a given composition;
* determine the phase decomposition of a phase that does not lie on the convex hull.

## Installation

* Navigate to the main directory.
* Make sure that the `setup.py` file is present.
* Run: `pip install .` (make sure to run the command in the desired Python environment).

## Structure

The code is implemented as a Python class named `cvhull`, which provides several attributes and methods. The class is documented, so its usage and that of its methods can be accessed using Python's built-in `help()` function.

The directory also contains a tutorial Jupyter notebook demonstrating the main features of the package (`tutorial.pynb`).

## Dependencies

* `numpy` (inclued in installation)
*  `pandas` (optional, but used in the tutorial)

## `cvhull` class

The `cvhull` class takes a set of N phases containing c elements as input and determines which phases are stable and unstable.

### Parameters

* `name` str[N] : list of names of the input phases
* `composition` float[N,c] : array containing the compositions of the phases
* `energy` float[N] : array containing the energies of the phases
* `species` str[c] : list of chemical species (default: A, B, C, ...)
* `accuracy` int : number of decimal places used for the compositions (default: 4)


*All phases must have unique names.*
*The code requires all reference phases (pure elements) to be provided. They must therefore be included in the input, even if their energy is equal to 0.*

### Attributes

* `cvhull.name` str[N] : list of names of the input phases
* `cvhull.composition` float[N,c] : array containing the phase compositions
* `cvhull.energy` float[N] : array containing the phase energies
* `cvhull.npoints` int : number of phases
* `cvhull.dim` int : dimensionality of the phase space
* `cvhull.accuracy` int : number of decimal places used for the compositions

* `cvhull.stable` int[Ns] : index of the stable phases
* `cvhull.instable` int[Ni] : index of the unstable phases

For each of the nf facets of the convex hull:                                              

* `cvhull.nfacets` int : number of facets
* `cvhull.hplan` float[nf,c+1] : parameters (p_i) defining the hyperplane equation (e.g. in 2D: $p_1x+p_2y+p_3=0$)
* `cvhull.vertix float[nf,c+1] : Vertices of the facet

### Methods

* `cvhull.energy_hull(compo)`
  Return the energy of the convex hull at the specified composition `compo`

* `cvhull.energy_compound(compound)`
  Return the energy of the phase `compound`

* `cvhull.print_stables()`
  Print the list of stable phases

* `cvhull.equilibrium(compo)`
 Display the phases equilibrium at composition `compo`

* `cvhull.distance_hull(name)`
  Display the deviation from the convex hull and the decompotion phases

* `cvhull.info_compound(name)`
 Display all info about the compound `name`

* `cvhull.table_dhull(name)`
  Return a dictionory of stability for all compound
