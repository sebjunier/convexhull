# -*- coding: utf-8 -*-
"""
Calcul d'enveloppe convexe, au sens thermodynamiques, à n dimensions.
"""

#Infos code
VERSION    = 3.1 # Date : 30/05/2023
AUTHOR     = 'S. Junier'

#Module
from scipy.spatial import ConvexHull as cvh
import numpy as np
import sys
from .fonction import *

class cvhull :
    """
Déterminer la stabilité des phases par calcul d'enveloppe convexe. Elle prend en entrée une nombre N de phases de c éléments

Paramètres :
------------
nom ([str], len=N) : Nom des phases
composition (np.array([N,c], flaot) : Compositions des phases
energie (np.array([N], float) : Energie des phases
accuracy (int, default=4) : Nombre de chiffres après la virgule utilisé pour les compositions

Les entrées np.array peut être mise sous forme de liste équivalente.

Attributs :
-----------
nom ([str], len=N) : Nom des phases
composition (np.array([N,c],float)) : Compsitions des phases
energie (np.array([N],float)) : Energie des phases
points (np.array([N,c],float)) : Points dans le format quickhull [composition sauf une + energie]
npoints (int) : Nombre de phases
dim (int) : Dimension des phases
accuracy (int) : Nombre de chiffres après la virgule utilisé pour les compositions

stable (np.array([N],int)) : Index des phases stables
instable (np.array([N],int)) : Index des phases  instables

Pour les nf facettes de l'enveloppe convexe :
equations (np.array([nf,c+1],float)) : Paramètres (pi) de l'équation de l'hyper-plan (\sum_{i=0}^{N-3}(p[i] * x[i]) + p[N-2]*y + p[N-1] = 0 )
sommets (np.array([nf,c+1],float)) : Sommets de la facette

"""
    def __init__(self, nom, composition, energie, accuracy=4) :
        #========================================= Verification des données d'entrée ==============================================================
        #Tableau des compositions
        try :
            np.array(composition, dtype='float64')
        except ValueError :
            print("Erreur dans les paramètres d'entrée :")
            print("La liste des 'composition' doit être de dimension fixe")
            sys.exit()
        for compo in composition :
            somme = 0
            for x in compo :
                somme += x
            if round(somme, accuracy) != 1 :
                print("Erreur dans les paramètres d'entrée :")
                print("La somme de chacun  des éléments de 'composition' doit être égale à 1")
                sys.exit()
                
        #Nombre de points
        if len(nom) <= len(composition[0]) :
            print("Erreur : Il faut au minimum un nombre de points superieur à la dimensions des points")
            print("Ici il faut donc au minimum {} points alors que vous en avez mit {}".format(len(composition[0])+1, len(nom)))
            sys.exit()
        
        #Tableau des noms
        tmp = []
        if type(nom) != list :
            print("Erreur dans les paramètres d'entrée :")
            print("la variable 'nom' doit être une liste")
            sys.exit()
        for x in nom :
            if type(x) != str :
                print("Erreur dans les paramètres d'entrée :")
                print("'nom' doit être une liste de str")
                sys.exit()
            if not x in tmp :
                tmp.append(x)
            else :
                print("Erreur dans les paramètres d'entrée :")
                print("la liste 'nom' doit avoir des éléments différents")
                print(x, 'est deux fois dans la liste')
                sys.exit()
        del tmp
        
        #Tableau des énergies
        try :
            np.array(energie, dtype='float64')
        except :
            print("Erreur dans les paramètres d'entrée :")
            print("'energie' doit être une liste de nombre")
            sys.exit()
            
        #Taille des tableaux d'entrée
        if not len(nom) == len(composition) == len(energie) :
            print("Erreur dans les paramètres d'entrée :")
            print("Les variables 'nom', 'composition' et 'energie' doivent avoir la même longeur")
            print("ici : ", len(nom) ,len(composition), len(energie))
            sys.exit()
            
        #===============================================================================Attributs===========================================================
        
        self.nom         = nom                                    # nom des composés
        self.composition = np.array(composition, dtype='float64') # compositions
        self.energie     = np.array(energie, dtype='float64')     # energies
        self.accuracy    = accuracy                               # Nombre de chiffre apres la virgule pour les compositions
        self.npoints     = len(composition)                       # nombre de points
        self.dim         = len(composition[0])                    # Dimension des points
        
        #On cherche l'energie des élements purs
        Epur = [None for _ in range(self.dim)]             # Energie du l'élément pur sur l'axe i
        index_pur  = [None for _ in range(self.dim)]       # indice de l'éléments pur sur l'axe i (dans le tableau self.nom)
        check = np.array([False for _ in range(self.dim)]) # True si on a deja un point sur l'axe i False sinon
        
        for i in range(self.npoints) :
            for j in range(self.dim) :
                if self.composition[i,j] == 1 : # si la composition de l'élément j, alors le composé est un élément pur de type j
                    if not check[j] : # si il n'y a acun point sur l'axe j, on ajoute le point i
                        check[j] = True
                        Epur[j] = self.energie[i]
                        index_pur[j] = i

                    elif self.energie[i] < Epur[j] : # sinon on l'ajoute seulement si sont énergie est inférieur
                        Epur[j] = self.energie[i]
                        index_pur[j] = i

        if not check.all() : # On vérfie que tout les élements purs sont dans la base
            sys.exit("Manque certains paramètres : Il faut mettre les éléments purs")
         
        # Pour le calcul de l'enveloppe convexe, on ne prend que les points d'énergie plus basse que la combinaison des élements purs
        points  = [] # Points au format quickhull [composition sauf une + energie]
        indices = [] # indice des points dans le tableau self.nom
        for i in range(self.npoints) :
            if self.energie[i] <= hyperplan_Epur(self.composition[i], Epur) :
                points.append(format_quickhull(composition[i], energie[i]))
                indices.append(i)
                
        #On vérifie que les elements purs ont bien été pris
        for i in index_pur :
            if not i in indices :
                points.append(format_quickhull(composition[i], energie[i]))
                indices.append(i)

        #Si les seuls points stables sont les éléments purs ont écrit les tableaux equations et sommets à la main
        if len(points) == self.dim : 
            print('-----------------------------------')
            print('Seul les composés purs sont stables')
            print('-----------------------------------')
            eq = []
            #paramètres p de l'equation d'un hyper plan passant par les ndim points Epur : \sum_{i=0}^{N-3}(p[i] * x[i]) + p[N-2]*y + p[N-1] = 0
            for i in range(self.dim-1) :
                eq.append(Epur[-1]-Epur[i])
            eq.append(1)
            eq.append(-Epur[-1])
            
            self.equations = [eq]
            self.sommets   = [index_pur]
            self.stable    = np.array(index_pur,dtype='int32')
            
            instable = []
            for i in range(self.npoints) :
                if not i in self.stable :
                    instable.append(i)
            self.instable = np.array(instable, dtype='int32')
        
        #Si il y a plus de composé stable
        else:
            qh = cvh(points)
            #On cherche quels sont les points stables
            stable   = []
            instable = []
            
            for i in qh.vertices :
                stable.append(indices[i]) #Les indices de qh ne sont pas bons
                                          #L'indice correct des points est dans le tableau 'indices'
            self.stable = np.array(stable, dtype='int32') # Liste des index des composés stables
            for i in range(self.npoints) :
                if not i in self.stable :
                    instable.append(i)
            self.instable = np.array(instable, dtype='int32') # Liste des index des composés pas stable
            
        #Selection des sommets et des facttes
            equation    = [] # Parametre p de l'equation de chaque facette  : \sum_{i=0}^{N-3}(p[i] * x[i]) + p[N-2]*y + p[N-1] = 0 
            sommets     = [] # indices des sommets de chaque facette
            ener_centre = [] # energie de chaque facette au point équimolaire
            x           = 1/self.dim
            x           = [ x for _ in range(self.dim-1) ] # point equimolaire
            
            #On prens toutes les facettes de l'enveloppe sauf celles orthogonale à l'énergie
            for i in range(len(qh.equations)) :
                if round(qh.equations[i][-2], 10) != 0 : # plans orthogonal à l'energie
                    equation.append(qh.equations[i])
                    sommets.append([ indices[j] for j in qh.simplices[i] ])
                    ener_centre.append(f_eq(x, qh.equations[i]))
            
            #On enleve la facette de plus haute énergie (partie haute de l'enveloppe)
            index_pop = ener_centre.index(max(ener_centre))
            equation.pop(index_pop)
            sommets.pop(index_pop)
            self.equations = equation
            self.sommets   = sommets

    #==============================================================================Méthode d'instance===========================================================  
    def energie_hull(self, compo) :
        """Donne l'énergie de l'enveloppe convexe à la composition 'compo'

        Args:
            compo ([float], len=self.dim): composition

        Returns:
            float : energie de l'enveloppe convexe à la composition 'compo'
            int   : indice de l'hyper plan correspondant
        """
        #Equation d'un hyperplan à partir des paramètres p (exemple : p1x + p2y + p3 = 0)
        #Verification de la composition
        somme = 0
        for x in compo :
            somme += x
        if round(somme, self.accuracy) != 1 :
            print("Erreur dans les paramètres d'entrée :")
            print("La somme des éléments de 'composition' doit être égale à 1")
            sys.exit()
        #Méthode
        E = []
        for param in self.equations :
            E.append(f_eq(compo[:-1], param))
        ieq = E.index(max(E)) # indice le l'hyper-plan corespondant
        return max(E), ieq 
    
    def energie_compose(self, nom):
        """
        Méthode qui donne l'énergie du composé 'nom'
        """
        index = self.nom.index(nom)
        return self.energie[index]
    
    def phases_stables(self, compo, accuracy=None) :
        """
        Affiche les phases stables à la composition 'compo'
         
         Args:
            compo ([float], len=self.dim): composition
            accuracy (int) : nombre de chiffres après la virgule pour les compositions
        """
        if accuracy == None :
            accuracy = self.accuracy
            
        ehull,ieq = self.energie_hull(compo) # On extrait l'indice de l'hyper plan correspondants aux phases stables
        
        decomposition_points = self.sommets[ieq] # indices des phases stables
        
        nom_dcp = [ self.nom[i] for i in decomposition_points ] # nom des phases stables
        
        # on veux résoudre la matrice qui donne la proportion de chaques phases de 'decomposition_points' de tel sorte que la composition total soit égale à 'compo' 
        matrix  = np.ones([self.dim,self.dim+1]) # Matrice remplie de 1
         
        for j,pt in enumerate(decomposition_points) : # chaque colonne corespond à une phases stables
            for i in range(self.dim-1) :              # chaque ligne corespond à une dimension de la composition
                matrix[i,j] = self.composition[pt][i]
        
        for i in range(self.dim-1) :
            matrix[i,-1] = compo[i] # la dernière colonne est la compostion 'compo'
        
        #On laisse que des 1 à la dernière ligne car la somme de la proportion de chaque phases doit être égale à 1
        
        result = resol(matrix) # matrice après résolution
        rtest  = np.array([[ round(result[i,j],10) for i in range(self.dim) ] for j in range(self.dim) ])
        if (rtest == identite(self.dim)).all() : # on vérifie que la dernière colonne représente bien la proportion de chaque phase 
            proportion = result[:,-1] # proportion de chaques phases
            
            fmt2 = "{:."+str(accuracy)+"f}({:}) + "
            print("Les phases stables sont :")
            txt = ""
            for i in range(self.dim) :
                if proportion[i] != 0 :
                    txt += fmt2.format(proportion[i],nom_dcp[i])
            txt = txt[:-3]
            print(txt) # On affiche la proportion de chaques phases
            #return(result)
            
        else :
            print(rtest)
            print(identite(self.dim))
            sys.exit('Erreur dans la resolution des équations')

    def distance_hull(self, nom, accuracy=None, sortie=True):
        """
        Méthode qui donne l'écart en énergie par rapport à l'enveloppe convexe du composé 'nom'
   
        Args:
            nom (str) : nom du composé
            accuracy (int) : nombre de chiffres après la virgule pour les compositions
            sortie (bool) : si True Affiche la décomposition de phases
            
        Returns:
            float : Différence d'énergie entre le composé 'nom' et l'enveloppe convexe
            bool : Stabilité du composé 'nom'
        """
        if accuracy == None :
            accuracy = self.accuracy
        
        index       = self.nom.index(nom)
        compo       = self.composition[index]
        ehull,ieq   = self.energie_hull(compo) # energie de l'enveloppe convexe à la compostion de 'nom'
        
        decomposition_points = self.sommets[ieq]
        nom_dcp     = [self.nom[i] for i in decomposition_points]
        
        ecompo      = self.energie[index] #enrgie du composé
        de          = ecompo-ehull # différence d'énergie
        
        if round(de, accuracy) == 0 :
            if sortie :
                print("Le composé",nom," est stable")
            ifstable = True
        else :   # Si le composé n'est pas stable on cherche la décomposition de phases (meme calcul que self.phases_stables)
            ifstable = False
            if sortie :
                matrix = np.ones([self.dim,self.dim+1])
                for j,pt in enumerate(decomposition_points) :
                    for i in range(self.dim-1) :
                        matrix[i,j] = self.composition[pt][i]
                for i in range(self.dim-1) :
                    matrix[i,-1] = compo[i]
                result = resol(matrix)
                if (result[:,:-1] == identite(self.dim)).all() :
                    proportion = result[:,-1]
                    fmt        = "{:."+str(accuracy)+"f}"
                    fmt2       = "{:."+str(accuracy)+"f}({:}) + "
                    
                    print("Le composé",nom," n'est pas stable")
                    print("--------------------------------------------")
                    print("\u0394H = "+fmt.format(de))
                    print("Il se décompose en :")
                    txt=""

                    for i in range(self.dim) :
                        if proportion[i] != 0 :
                            txt += fmt2.format(proportion[i],nom_dcp[i])
                    txt = txt[:-3]
                    print(txt)

                else :
                    print(result[:,:-1])
                    print(identite(self.dim))
                    sys.exit('Erreur dans la resolution des équations')
                    
        return de,ifstable
    
    def tableau_complet(self, elements_purs, save_csv=False) :
        """Extrait dans un taleau la différence par rapport à l'enveloppe convexe pour chaque composé

        Args:
            elements_purs ([str], len=dim): Liste des éléments du tableau composition
            save_csv (bool, optional): Si True, extrait le tableau dans un fichier data.csv (Defaults to False)

        Returns:
           list : clés du tableau
           lsit : valeurs du tableau 
        """
        keys_tab = ['Nom']
        
        if len(elements_purs) != self.dim :
            sys.exit('elements_purs doit etre une liste de la même taille que la dimension du problème')
        for elem in elements_purs :
            if type(elem) == str:
                keys_tab.append('x('+elem+')')
            else :
                sys.exit('elements_purs doit etre une liste de chaine de caractère')
        keys_tab.append('Stabilité')
        keys_tab.append("Ecart à l'enveloppe")
        
        N_keys = len(keys_tab)
        
        out_tab = []
        for i in range(self.npoints) :
            out_tab.append([])
            out_tab[-1].append(self.nom[i])
            
            for j in range(self.dim) :
                out_tab[-1].append(self.composition[i,j])
            diff_hull, ifstable = self.distance_hull(self.nom[i],sortie=False)
            if ifstable :
                out_tab[-1].append('oui')    
                out_tab[-1].append(0)
            else :
                out_tab[-1].append('non')    
                out_tab[-1].append(diff_hull)
        
        if save_csv :
            with open('data.csv','w') as file :
                line = ''
                for x in keys_tab :
                    line += x + ','
                file.write(line[:-1]+'\n')
                for i in range(self.npoints) :
                    line = ''
                    for x in out_tab[i]:
                        line += str(x) + ','
                    file.write(line[:-1]+'\n')
        
        return keys_tab, out_tab 
        
            
        
