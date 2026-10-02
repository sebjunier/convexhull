# -*- coding: utf-8 -*-

#Module
import numpy as np
from .function import *

class cvhull :
	"""
Class to calculate convexhull for materials data.


Parameters :
------------
* name (list, len=N) : Names of the compounds
* composition (np.array([N,c], float)) : Compositions of the compounds (or number of atoms)
* energy (np.array([N], float)) : Energy of the phases
* species (list(str), len=c) : List of chemical species (default: A, B, C, ...)
* accuracy (int) : Number of significant digits used for the compositions (default: 4)

N: number of compound
c: number of species

Attributs :
-----------
* name (list, len=N) : Name of compounds
* composition (np.array([N,c], float) : Composition of compounds
* energy (np.array([N], float) : Energy of the phases
* species (list, len=c) : list of chemical species (default: A,B,C,...)
* npoints (int) : Compounds number
* dim (int) : Species number
* accuracy (int) : Significatif numbers of compositions

* stable (np.array([N],int)) : Index of compounds on the convex hull
* instable (np.array([N],int)) : Index of compounds out of the convex hull

For each of the nf facets of the convex hull:
* nfacets (int) : number of facets
* hplan (np.array([nf,c+1], float)) : Parameters (p_i) defining the hyperplane equation
* vertix (np.array([nf,c+1], float)) : Vertices of the facet

"""
	def __init__(self, name, composition, energy, species=None, accuracy=4) :
		name, composition, energy = checkParameters(name, composition, energy)
		self.name = name
		self.composition = composition
		self.energy = energy
		self.accuracy = accuracy
		self.npoints  = composition.shape[0]
		self.dim      = composition.shape[1]
		if type(species) ==  type(None) :
			self.species = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:self.dim])
		else :
			self.species = species

		vertix = np.zeros((0,self.dim),dtype='int')
		hPlan  = np.zeros((0,self.dim+1),dtype='int')
		# Found all facets
		vertix_ref = getRef(composition, energy, self.dim,self.npoints)
		vertix_ref.sort()
		i = 0
		done = 0
		vertix, hPlan = addFacet(vertix_ref, composition, energy, vertix, hPlan)
		facetTested = vertix.copy() # list of already tested facets

		while done < 1 :
			duplhp = np.dot(np.ones((self.npoints,1)), hPlan[i].reshape((1,-1)))
			distance = energy - np.array(list(map(fHP, composition, duplhp)))
			index = np.argmin(distance)

			if distance[index] > -1e-12 :
				i += 1 # This facets is good, go to next
			else :
				for j in range(vertix.shape[1]) :
					new = np.concatenate((vertix[i,:j], vertix[i,j+1:], [index]))
					new.sort()
					if (new == facetTested).all(1).any() : # already tested
						continue
					facetTested = np.concatenate([facetTested,[new]])
					vertix, hPlan = addFacet(new, composition, energy, vertix, hPlan)
					
				vertix = np.delete(vertix, i, axis=0)
				hPlan = np.delete(hPlan, i, axis=0)
			done = i/len(vertix)

		self.hplan = hPlan
		self.vertix = vertix
		self.nfacets = vertix.shape[0]
		self.stable = np.array(list(set(vertix.reshape(-1))))
		self.stable.sort()
		self.instable = np.array([ x for x in range(self.npoints) if x not in self.stable])
    
#================================== Instance methods ===========================================  
	def energie_hull(self, compo, returnFacet=False) :
		"""Return energy of convexhull at composition 'compo' """
		compo = np.array(compo, dtype='float')
		compo = compo/compo.sum()

		duplcompo = np.dot(np.ones((self.nfacets,1)), compo.reshape((1,-1)))
		energyList = np.array(list(map(fHP, duplcompo, self.hplan)))
		index = np.argmax(energyList)
		if returnFacet :
			return energyList[index], index
		else :
			return energyList[index]
    
	def energy_compound(self, compound):
		"""Return energy of copound `compound`"""
		if compound not in self.name :
			raise ValueError(f"{compound} is not in the data abse")
		index = self.name.index(compound)
		return self.energy[index]

	def print_stables(self) :
        """ Print the list of stable phases."""
		print(f'The database consist of {self.npoints} compounds with {len(self.stable)} stable')
		print()
		print('List of stable compounds:')
		print('-------------------------')

		i = len(str(self.npoints-1))
		n = len(str( max(self.name, key=lambda x:len(str(x))) ))
		c = self.accuracy
		e = len(str(int(max(np.fabs(self.energy)))))+1
		e = max(3,e)
		fmt_title = f' [:>{i}]  ' + f' [:{n}]  ' + self.dim*f' [:>{c+3}]'     + f'   [:>{e+3}]'
		fmt       = f' [:>{i}]  ' + f' [:{n}]  ' + self.dim*f' [:{c+3}.{c}f]' + f'   [:{e+3}.2f]'
		line = (i+3+ n+3 + self.dim*(c+4) + e+6)*'-'

		fmt_title = fmt_title.replace('[','{').replace(']','}')
		fmt = fmt.replace('[','{').replace(']','}')
		print(line)
		print(fmt_title.format('i', 'Name', *self.species, 'Energy'))
		print(line)
		for i in self.stable :
			print(fmt.format(i, self.name[i], *self.composition[i,:], self.energy[i]))
		print(line)

	def equilibrium(self, compo) :
		"""Display the phases equilibrium at compostion `compo` """
		compo = np.array(compo, dtype='float')
		compo = compo/compo.sum()

		ehull, ieq = self.energie_hull(compo, returnFacet=True) # Hull energy and facet index
		dp_index = self.vertix[ieq] # stable phases (Decomposion Points)
		dp_name = [ self.name[i] for i in dp_index ]
		# Resolve matrix of composition
		matrix  = np.eye(self.dim) #square matrix
		for j,i in enumerate(dp_index) :
			matrix[:,j] = self.composition[i,:]

		result = np.linalg.solve(matrix, compo)
		order = sorted(range(self.dim), key=lambda x:result[x],reverse=True)
		fmt = f'[:.{self.accuracy}]([:])'
		fmt = fmt.replace('[','{').replace(']','}')
		print("Decomposition:", end=" ")
		step=""
		for i in order :
			if result[i] > 1e-12 :
				print(step, fmt.format(result[i], dp_name[i]), end=" ")
				step='+'
		print()

	
	def distance_hull(self, name, verbose=True):
		"""Display the deviation from the convex hull and the decompotion phases"""
		index = self.name.index(name)
		compo = self.composition[index,:]
		ehull = self.energie_hull(compo)
		if index in self.stable :
			dhull = 0
		else:
			dhull = self.energy[index]-ehull

		if verbose :
			if dhull :
				fmt = f'[:.{self.accuracy}]'
				fmt = fmt.replace('[','{').replace(']','}')

				print(f"The compound {name} is not stable")
				print("--------------------------------------------")
				print("\u0394H = "+fmt.format(dhull))
				self.equilibrium(compo)
			else :
				print(f"The compound {name} is stable!")
                    
		return dhull

	def info_compond(self, name) :
        """Display all info about the compound `name` """
		print(f'Info about compound : {name}')
		print('------------------------------------')
		print(f'E = {self.energy_compound(name)}')
		self.distance_hull(name)


	def table_dhull(self) :
		""" Return a dictionory of stability for all compound"""
		keys_tab = ['Name']+list(self.species)+["Stability", "Deviation from hull"]
		out_tab = {}
		for i in keys_tab :
			out_tab[i] = []

		for i in range(self.npoints) :
			out_tab['Name'].append(self.name[i])

			for j in range(self.dim) :
				index = self.species[j]
				out_tab[index].append(self.composition[i,j])

			dhull = self.distance_hull(self.name[i], verbose=False)

			out_tab["Stability"].append(not dhull) # If stable True, esle False
			out_tab["Deviation from hull"].append(dhull)

		return out_tab
