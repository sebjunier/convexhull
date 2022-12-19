#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Module
import numpy as np
import sys

#----------------------------------------------------------------------------------Fonctions
#Equation d'un hyperplan à partir des paramètres p (exemple : p1x + p2y + p3 = 0)
def f_eq(x,p):
    if len(x)!=len(p)-2:
        sys.exit('le nombre de variables x doit etre égale au nombre de paramètre p - 2')
    e=-p[-1]/p[-2]
    for i in range(len(x)):
        e-=x[i]*p[i]/p[-2]
    return e
#Fonction de l'hyper-plan des elemetns pures (à partir des E de chacun des éléments purs)
def fp(xs,es):
    if len(xs)!=len(es)-1:
        sys.exit('le nombre de variables xs doit etre égale au nombre de points es -1')
    e0=es[-1]
    e=es[-1]
    for i  in range(len(xs)) :
        e+=(es[i]-e0)*xs[i]
    return e
#Resoud le systeme y ; la derniere colones de y corespond au terme à droite du signe égal
def resol(y,acc=10):
    x=np.array(y,float)
    
    variables=[]
    N=len(x)
    l_ligne=list(range(N))
    
    
    for C in range(N):
        flag=True
        for L in l_ligne :
            if x[L,C]!=0:
                flag=False
                l_ligne.remove(L)
                variables.append(L)
                for lin in range(N):
                    if lin!=L:
                        fact=x[lin,C]/x[L,C]
                        for col in range(N+1):
                            x[lin,col]-=fact*x[L,col]
        
                break
        
        if flag :
            variables.append(None)
    
    flag2=False        
    
    for lin in l_ligne:
        if x[lin,-1]!=0:
            flag2=True
    
    if flag2:
        print('Certaines equation se contredise')
    else:
        x_new=[]
        for C,lin in enumerate(variables):
            if lin!=None:
                fact=1/x[lin,C]
                for col in range(N+1):
                    x[lin,col]*=fact
                x_new.append(x[lin])
    
    x=np.array(x_new)
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            x[i,j]=round(x[i,j],acc)
    return x
#Ecrit la marice identité de dimention x
def identite(x):
    id=np.zeros([x,x])
    for i in range(x):
        id[i,i]=1
    return id
