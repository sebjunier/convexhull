# -*- coding: utf-8 -*-
"""
Principals functions of package convexhull
"""
#Module
import numpy as np

def getHP(matrix) : # Get Hyper-plan parameters from matrix
	m = matrix.shape[1] - 1
	return np.linalg.solve(matrix, np.concatenate([np.zeros((m)),[1]]))

def fHP(compo, hp) : #Return values of convex hull at composition `compo`
	return -(np.dot(compo[:-1], hp[:-2]) + hp[-1])/hp[-2]

def addFacet(new, composition, energy, vertix, hPlan, tol=1e-12) :
	m = composition.shape[1]
	matrix = np.zeros((m+1,m+1))
	points = np.zeros((m,m))
	j = 0
	for i in new :
		points[j,:] = composition[i,:]
		matrix[j,:-2] = composition[i,:-1]
		matrix[j,-2] = energy[i]
		matrix[j,-1] = 1
		j+=1
	matrix[-1,-2] = 1

	rank = np.linalg.matrix_rank(points)

	if rank == m :
		try :
			hp = getHP(matrix)
		except :
			print("Error")
			return vertix, hPlan
		vertix = np.concatenate([vertix, [new]])
		hPlan = np.concatenate([hPlan, [hp]])
	return vertix, hPlan

def getRef(composition, energy, dim, N) :
	indexList = np.arange(dim)

	equilibrium = np.zeros((dim),dtype='int')

# Look for the more stable reference
	ref = np.array([None]*dim)
	for i in range(N) :
		x = composition[i,:]
		if (x==1).any() :
			j = int((indexList*x).sum())
			if ref[j] == None or ref[j]>energy[i] :
				ref[j] = energy[i]
				equilibrium[j] = i
	if (ref==None).any() : # Check il there is a compound for each reference
		raise ValueError('All reference (purs elements) must be in database')
	return equilibrium

def checkParameters(name, composition, energy) :
	try :
		composition = np.array(composition, dtype='float64')
	except ValueError :
		raise ValueError("Composition must be a matrix of N,m of float")
	composition = composition/composition.sum(axis=1, keepdims=True)

	try :
		if len(set(name)) != len(name) :
			raise ValueError("name must have unique elements")
	except TypeError :
		raise TypeError("name must be iterable")

	try :
		energy = np.array(energy, dtype='float64')
	except :
		raise TypeError("energie must be a list of float")

	if not len(name) == len(composition) == len(energy) :
		raise ("The arguments `name`, `composition` and `energy` must have the same size")

	return list(name), composition, energy

