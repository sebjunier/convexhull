#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Infos code
VERSION=2.3
NAME_CODE='cvh'
DEPENDANCE='cvh_fonc.py'

#Module
from scipy.spatial import ConvexHull as cvh
import numpy as np
import sys
from cvh_fonc import *

#Classe
class cvhull :
    """
Classe pour déterminer la stabilité des composé par calcul d'enveloppe convexe
----------------------------------------------------------------------------------------------------
La classe prend en entrée une nombre N de composé de c éléments et caculs les composé stable et
instable

Paramètres:
-----------
nom | [str(N)]                : liste des nom des composés d'entrée
composition | [[float(c)](N)] : liste des compsitions des composés
energie | [flaot(N)]          : liste des énergie des composés

decimale | int (default=4)    : nombre de chiffres après la virgule utlisés pour les compositions

Atributs:
---------
nom | [str(N)] : Liste des nom des composés d'entrée

composition | np.array([N,c],float) : Liste des compsitions des composés
energie     | np.array([N],float)   : Liste des énergie des composés
points      | np.array([N,c],float) : Liste de points dans le format quickhull 
                                      [composition sauf une + energie]
npoints  | int : Nombre de points
dim      | int : Dimension des composés
decimale | int : Nombre de chifres après la vigures considéré

stable   | np.array([N],int) : Liste des index des composés stables
instable | np.array([N],int) : Liste des index des composés pas stable

equations | np.array([n,c+1],float) : Liste des paramètres (pi) des équations des plans de 
                                      l'enveloppe convexe (ex à 2D: p1x + p2y + p3 = 0)
sommets   | np.array([n,c+1],float) : Liste des sommets de l'hyper-plan
Méthodes:
---------
energie_hull(compo)
energie_compose(nom)
distance_hull(nom)
"""
    def __init__(self,nom,composition,energie,decimale=4) :
        #-----------------------------------------------------------Verification des donnée d'entrée
        #Taille des tableau d'entré
        if not len(nom)==len(composition)==len(energie):
            print("Erreur dans les paramètres d'entrée :")
            print("Les variables 'nom', 'composition' et 'energie' doivent avoir la même longeur")
            print("ici : ",len(nom),len(composition),len(energie))
            sys.exit()
        #Tableau des compositions
        try :
            np.array(composition,dtype='float64')
        except ValueError :
            print("Erreur dans les paramètres d'entrée :")
            print("La liste des 'composition' doit être de dimension fixe")
            sys.exit()
        for compo in composition:
            somme=0
            for x in compo:
                somme+=x
            if round(somme,decimale)!=1:
                print("Erreur dans les paramètres d'entrée :")
                print("La somme des éléments de 'composition' doit être égale à 1")
                sys.exit()
        #Nombre de points
        if len(nom)<=len(composition[0]):
            print("Erreur : Il faut au minimum un nombre de points superieur à la dimensions des points")
            print("Ici il faut donc au minimum {} points alors que vous en avez mit ".format(len(composition[0])+1))
            sys.exit()
        #Tableau des noms
        tmp=[]
        if type(nom)!=list:
            print("Erreur dans les paramètres d'entrée :")
            print("la variable 'nom' doit être une liste")
            sys.exit()
        for x in nom:
            if type(x)!=str:
                print("Erreur dans les paramètres d'entrée :")
                print("'nom' doit être une liste de str")
                sys.exit()
            if not x in tmp:
                tmp.append(x)
            else :
                print("Erreur dans les paramètres d'entrée :")
                print("la liste 'nom' doit avoir des éléments différents")
                print(x,'est deux fois dans la liste')
                sys.exit()
        del tmp
        #Tableau des énergies
        try :
            np.array(energie,dtype='float64')
        except :
            print("Erreur dans les paramètres d'entrée :")
            print("'energie' doit être une liste de nombre")
            sys.exit()
        #----------------------------------------------------------------------------------Attributs
        self.nom=nom # nom des composés
        self.composition=np.array(composition,dtype='float64') # composition
        self.energie=np.array(energie,dtype='float64') # energies
        self.decimale=decimale #Nombre de chiffre apres la virgule pour les composition
        self.npoints=len(composition) #nombre de points
        self.dim=len(composition[0]) # Dimension des points
        
        #On cherche l'energie des élements pur
        E0s=[True for _ in range(self.dim)] # Energie du l'élément pur sur l'axe i
        check=np.array([False for _ in range(self.dim)]) # True si on a deja un points sur l'axe i
        iE0s=[None for _ in range(self.dim)] # indice de l'éléments pur sur l'axe i
        for i in range(self.npoints):
            for j in range(self.dim):
                if self.composition[i,j]==1:
                    if not check[j]:
                        check[j]=True
                        E0s[j]=self.energie[i]
                        iE0s[j]=i
                    elif self.energie[i]<E0s[j]:
                        E0s[j]=self.energie[i]
                        iE0s[j]=i

        if not check.all(): 
            sys.exit("Maque certains paramètres : Il faut mettre les éléments purs")
         
        #Pour le calcul de l'enveloppe convexe, on ne prend que les points
        #inférieur à la combinaison des élements purs
        points=[] #Points au format quickhull [composition sauf une + energie]
        indices=[] #indice des points dans le tableau self.nom
        for i in range(self.npoints):
            if self.energie[i]<=fp(self.composition[i][:-1],E0s):
                points.append([x for x in composition[i]][:-1]+[energie[i]])
                indices.append(i)
        #On vérifie que les elements purs ont bien été pris
        for i in iE0s:
            if not i in indices:
                points.append([x for x in composition[i]][:-1]+[energie[i]])
                indices.append(i)

        #Si les seuls points stables sont les éléments purs
        if len(points)==self.dim: 
            print('-----------------------------------')
            print('Seul les composés purs sont stables')
            print('-----------------------------------')
            eq=[]
            #Equation d'un hyper plan passant par les ndim points E0s
            for i in range(self.dim-1):
                eq.append(E0s[-1]-E0s[i])
            
            eq.append(1)
            eq.append(-E0s[-1])
            self.equations=[eq]
            self.sommets=[iE0s]
            self.stable=np.array(iE0s,dtype='int32')
            instable=[]
            for i in range(self.npoints):
                if not i in self.stable:
                    instable.append(i)
            self.instable=np.array(instable,dtype='int32')
        
        #Si il y a plus de composé stable
        else:
            qh=cvh(points)
            #On cherche quels sont les points stables
            stable=[]
            instable=[]
            for i in qh.vertices :
                stable.append(indices[i]) #Les indices de qh ne sont pas bons
                                          #L'indice correct des points est dans le tableau 'indices'
            self.stable=np.array(stable,dtype='int32') # Liste des index des composés stables
            for i in range(self.npoints):
                if not i in self.stable:
                    instable.append(i)
            self.instable=np.array(instable,dtype='int32') # Liste des index des composés pas stable
            
        #Selection des sommets et des facttes
            equation=[]
            sommets=[]
            ener_centre=[]
            x=1/self.dim
            x=[x for _ in range(self.dim-1)]
            
            #On prens tout les facettes de l'envelope sauf celles ortho à E
            for i in range(len(qh.equations)):
                if round(qh.equations[i][-2],10)!=0: # plans orthogonal à l'energie
                    equation.append(qh.equations[i])
                    sommets.append([indices[j] for j in qh.simplices[i]])
                    ener_centre.append(f_eq(x,qh.equations[i]))
            #On enleve la facette de plus haute énergie (partie haute de l'enveloppe)
            index_pop=ener_centre.index(max(ener_centre))
            equation.pop(index_pop)
            sommets.pop(index_pop)
            self.equations=equation
            self.sommets=sommets

        
    def energie_hull(self,compo):
        """
        Méthode qui donne l'énergie de l'enveloppe convexe à la composition 'compo'
        """
        #Equation d'un hyperplan à partir des paramètres p (exemple : p1x + p2y + p3 = 0)
        #Verification de la composition
        somme=0
        for x in compo:
            somme+=x
        if round(somme,self.decimale)!=1:
            print("Erreur dans les paramètres d'entrée :")
            print("La somme des éléments de 'composition' doit être égale à 1")
            sys.exit()
        #Méthode
        E=[]
        for param in self.equations:
            E.append(f_eq(compo[:-1],param))
        ieq=E.index(max(E)) # indice le l'hyper-plan corespondant
        return max(E),ieq 
    
    def energie_compose(self,nom):
        """
        Méthode qui donne l'énergie de l'énergie du composé 'nom'
        """
        index=self.nom.index(nom)
        return self.energie[index]
    
    def phases_stables(self,compo,decimale=4):
        """
        Méthode qui donne les phases stable à la composition 'compo'
        les valeurs sont exprimé avec 'decimale' chifres apres la virgules (defaut :4)'
        """
        ehull,ieq=self.energie_hull(compo) #verif de la compo dans la méthode
        decomposition_points=self.sommets[ieq]
        nom_dcp=[self.nom[i] for i in decomposition_points]
        matrix=np.ones([self.dim,self.dim+1])
        for j,pt in enumerate(decomposition_points):
            for i in range(self.dim-1):
                matrix[i,j]=self.composition[pt][i]
        for i in range(self.dim-1):
            matrix[i,-1]=compo[i]
        result=resol(matrix)
        rtest=np.array([[round(result[i,j],10) for i in range(self.dim)] for j in range(self.dim)])
        if (rtest==identite(self.dim)).all() :
            proportion=result[:,-1]
            fmt2="{:."+str(decimale)+"f}{:} + "
            print("Les phases stables sont :")
            txt=""
            for i in range(self.dim):
                if proportion[i]!=0:
                    txt+=fmt2.format(proportion[i],nom_dcp[i])
            txt=txt[:-3]
            print(txt)
            return(txt)
            
        else :
            print(rtest)
            print(identite(self.dim))
            sys.exit('Erreur dans la resolution des équations')

    def distance_hull(self,nom,decimale=4,sortie=True):
        """
        Méthode qui donne l'écart en énergie par rapport à l'enveloppe covexe du composé 'nom'
        si il n'est pas stable elle donne la décomposition
        Les valeur sont donnée avec 'decimale' chiffres après la virgule 
        """
        index=self.nom.index(nom)
        compo=self.composition[index]
        ehull,ieq=self.energie_hull(compo)
        decomposition_points=self.sommets[ieq]
        nom_dcp=[self.nom[i] for i in decomposition_points]
        ecompo=self.energie[index]
        de=ecompo-ehull
        if round(de,decimale)==0:
            print('Ce composé est stable')
            ifstable=True
        else :
            ifstable=False
            if sortie:
                matrix=np.ones([self.dim,self.dim+1])
                for j,pt in enumerate(decomposition_points):
                    for i in range(self.dim-1):
                        matrix[i,j]=self.composition[pt][i]
                for i in range(self.dim-1):
                    matrix[i,-1]=compo[i]
                result=resol(matrix)
                if (result[:,:-1]==identite(self.dim)).all() :
                    proportion=result[:,-1]
                    fmt="{:."+str(decimale)+"f}"
                    fmt2="{:."+str(decimale)+"f}{:} + "
                    print("Le composé",nom," n'est pas stable")
                    print("--------------------------------------------")
                    print("\u0394H = "+fmt.format(de))
                    print("Il se décompose en :")
                    txt=""
                    for i in range(self.dim):
                        if proportion[i]!=0:
                            txt+=fmt2.format(proportion[i],nom_dcp[i])
                    txt=txt[:-3]
                    print(txt)
                else :
                    print(result[:,:-1])
                    print(identite(self.dim))
                    sys.exit('Erreur dans la resolution des équations')
        return de,ifstable
        
