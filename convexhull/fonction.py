# -*- coding: utf-8 -*-
"""
Fonctions principales du paquet
"""
#Module
import numpy as np
import sys

def format_quickhull(compo, energie) :
    """transforme la composition et l'energie d'un composé au format quickhull
    """
    compo_inde = [x for x in compo][:-1] # on transforme compo en tableau (si c'est un np.array) et on enlève le dernier élement
    return compo_inde + [energie]

def f_eq(compo, p) :
    """Equation d'un hyperplan à partir des paramètres p (equation pour N paramètres : \sum_{i=0}^{N-3}(p[i] * compo[i]) + p[N-2]*y + p[N-1] = 0)

    Args:
        compo ([float], len=N-2) : composition indépendantes
        p ([float], len=N)   : paramètre d'une equation de facette

    Returns:
        float : energie de l'hyper-plan à la composition compo
    """
    if len(compo) != len(p)-2 :
        sys.exit('le nombre de variables compo doit etre égale au nombre de paramètre p - 2')
    Ehp = -p[-1] / p[-2]
    for i in range(len(compo)) :
        Ehp -= compo[i] * p[i] / p[-2]
    return Ehp

def hyperplan_Epur(compo, Epur) :
    """ fonction de l'hyper plan à partir de l'energie de chaque élément pur

    Args:
        compo ([float], len=c): composition
        Epur ([float], len=c): energie de chaque élément pur

    Returns:
        float : energie de l'hyper-plan à la composition compo
    """
    if len(compo) != len(Epur) :
        sys.exit('le nombre de variables compo doit etre égale au nombre de points Epur')
    E0 = Epur[-1]
    Ehp = Epur[-1]
    for i  in range(len(compo) - 1) : 
        Ehp += (Epur[i] - E0) * compo[i]
    return Ehp

#Resoud le systeme y ; la derniere colones de y corespond au terme à droite du signe égal
def resol(y, acc=10) :
    x = np.array(y, float)
    
    variables = []
    N = len(x)
    l_ligne = list(range(N))
    
    
    for C in range(N) :
        flag = True
        for L in l_ligne :
            if x[L,C] != 0 :
                flag = False
                l_ligne.remove(L)
                variables.append(L)
                for lin in range(N) :
                    if lin != L :
                        fact = x[lin,C] / x[L,C]
                        for col in range(N+1) :
                            x[lin,col] -= fact * x[L,col]
        
                break
        
        if flag :
            variables.append(None)
    
    flag2 = False        
    
    for lin in l_ligne :
        if x[lin,-1] != 0 :
            flag2 = True
    
    if flag2 :
        print('Certaines equation se contredise')
    else:
        x_new = []
        for C, lin in enumerate(variables) :
            if lin != None :
                fact = 1 / x[lin,C]
                for col in range(N+1) :
                    x[lin,col] *= fact
                x_new.append(x[lin])
    
    x = np.array(x_new)
    for i in range(x.shape[0]) :
        for j in range(x.shape[1]) :
            x[i,j] = round(x[i,j], acc)
    return x
#Ecrit la marice identité de dimention x
def identite(x) :
    id=np.zeros([x,x])
    for i in range(x) :
        id[i,i] = 1
    return id