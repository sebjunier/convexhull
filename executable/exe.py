import sys
import numpy as np
from convexhull import cvhull

#Entrée
file_in = 'data_in.csv'
file_out = 'data_out'
element = ['Fe','Cu','Co','Ni']
key_name = 'id'
key_energy = 'E'

# fonction d'extraction
def extract_tab(tab, keys) :
    keys_list = tab[0,:]
    if type(keys) == str :
        index = np.where(keys_list == keys)[0][0]
        y = tab[1:,index]
    elif type(keys) == list :
        indexs = []
        for key in keys :
            indexs.append(np.where(keys_list == key)[0][0])
        y = tab[1:,indexs]
    else :
        sys.exit()
    return y

#Extraction donée

data = np.loadtxt(file_in,dtype=str,delimiter=',')

energy = extract_tab(data,key_energy)
energy = np.array(energy,'float64')

name = extract_tab(data,key_name)
name = [x.strip() for x in name]

compo = extract_tab(data,element)
compo = np.array(compo,'float64')
compo = [[y/x.sum() for y in x] for x in compo]
compo = np.array(compo,'float64')

#Calcul et sorite

hull = cvhull(name, compo, energy)
hull.tableau_complet(element,save_csv=file_out)