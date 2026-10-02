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

* `numpy`

## `cvhull` class

The `cvhull` class takes a set of N phases containing c elements as input and determines which phases are stable and unstable.

### Parameters

* `nom` str[N] : list of names of the input phases
* `composition` float[N,c] : array containing the compositions of the phases
* `energie` float[N] : array containing the energies of the phases
* `prec` int : number of decimal places used for the compositions (default: 4)
  *All phases must have unique names.*
  *The code requires all reference phases (pure elements) to be provided. They must therefore be included in the input, even if their energy is equal to 0.*

### Attributes

* `cvhull.nom` str[N] : list of names of the input phases

* `cvhull.composition` np.array([N,c],float) : array containing the phase compositions

* `cvhull.energie` np.array([N],float) : array containing the phase energies

* `cvhull.points` np.array([N,c],float) : array of points in the format required by Quickhull [independent compositions + energy]

* `cvhull.npoints` int : number of phases

* `cvhull.dim` int : dimensionality of the phase space

* `cvhull.prec` int : number of decimal places used for the compositions

* `cvhull.stable` np.array([N],int) : indices of the stable phases

* `cvhull.instable` np.array([N],int) : indices of the unstable phases

* `cvhull.equations` np.array([n,c+1],float) : parameters $p_i$ of the hyperplane equations defining the convex hull (e.g. in 2D: $p_1x+p_2y+p_3=0$)

* `cvhull.sommets` np.array([n,c+1],float) : vertices defining the hyperplanes

### Methods

* `cvhull.energie_hull(compo)`
  Returns the energy of the convex hull at the specified composition `compo` (float[c]).

* `cvhull.energie_compose(nom)`
  Returns the energy of the phase `nom` (str).

* `cvhull.phases_stables(compo,prec=4)`
  Returns the phase decomposition at the specified composition `compo` (float[c]). The resulting values are given with `prec` decimal places (default: 4).

* `cvhull.distance_hull(nom,prec=4,sortie=True)`
  Returns the energy above the convex hull for the phase `nom` (str). If the phase is unstable, the method also returns its phase decomposition. This behavior can be disabled by setting `sortie=False`.
