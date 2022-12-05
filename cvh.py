#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
        # On met les points au format quickhull [composition sauf une + energie] 
        points=[[x for x in composition[i]][:-1]+[energie[i]] for i in range(len(nom))]
        self.npoints=len(points) #nombre de points
        self.points=np.array(points) #liste de points dans le format quickhull
        self.dim=len(points[0]) # Dimension des points
        
        #On cherche l'energie des élements pur
        E0s=[True for _ in range(self.dim)]
        check=np.array([False for _ in range(self.dim)])
        iE0s=[None for _ in range(self.dim)]
        for i in range(self.npoints):
            for j in range(self.dim):
                if self.composition[i,j]==1:
                    check[j]=True
                    if E0s[j]:
                        E0s[j]=self.energie[i]
                        iE0s[j]=i
                    elif self.energie[i]<E0s[j]:
                        E0s[j]=self.energie[i]
                        iE0s[j]=i

        if not check.all(): 
            sys.exit("Maque certains paramètres : Il faut mettre les éléments purs")
        
        qh=cvh(points)
        #Les élements stable sont les éléments de l'enveloppe convexe inférieur à la 
        #combinaison des élements purs
        stable=[]
        for i in qh.vertices:
            if self.energie[i]<=fp(self.composition[i][:-1],E0s):
                stable.append(i)
        self.stable=np.array(stable,dtype='int32') # Liste des index des composés stables
        instable=[]
        for i in range(self.npoints):
            if not i in self.stable:
                instable.append(i)
        self.instable=np.array(instable,dtype='int32') # Liste des index des composés pas stable
         
        #On cherche la valeur du centre de l'hyper-plan des elements pur
        E0=np.array(E0s).mean() 
        equation=[]
        sommets=[]
        x=1/self.dim
        x=[x for _ in range(self.dim-1)]
        
        #On ne considére pas les plan orthogonal à l'energie
        #On prens seulement les hyper plan qui ont un centre plus bas que celui des élement
        for i in range(len(qh.equations)):
            if qh.equations[i][-2]!=0:
                if f_eq(x,qh.equations[i])<E0:
                    equation.append(qh.equations[i])
                    sommets.append(qh.simplices[i])
        if len(equation)!=0:
            self.equations=equation
            self.sommets=sommets
        else :
            eq=[]
            for i in range(1,self.dim):
                eq.append(E0s[-1]-E0s[i])
            
            eq.append(1)
            eq.append(-E0s[-1])
            self.equations=[eq]
            self.sommets=[iE0s]
        
    def energie_hull(self,compo):
        """
        Méthode qui donne l'énergie de l'enveloppe convexe à la composition 'compo'
        """
        #Equation d'un hyperplan à partir des paramètres p (exemple : p1x + p2y + p3 = 0)
        #Tests
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
        ehull,ieq=self.energie_hull(compo)
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

    def distance_hull(self,nom,decimale=4):
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
                print("Ce composé n'est pas stable")
                print("-----------------------------")
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
        
