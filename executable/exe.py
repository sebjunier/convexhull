#!/usr/bin/env python
import sys
import os
import tkinter as tk
import numpy as np
from convexhull import cvhull

def error(*txt) : # Ouvre une fenêtre d'erreur
    erreur = tk.Tk()
    erreur.geometry('700x200')
    erreur.config()
    erreur.title("Error")
    lentry = tk.Label(erreur, text="Erreur :", font=("arial",15))
    lentry.pack()
    for t in txt :
        lentry = tk.Label(erreur, text=t, font=("arial",12))
        lentry.pack()
    

def clic_bouton(): # Lance le calcul si tout les paramètre sont ok 
    global file_in, file_out, element, key_energy, key_name, flag
    
    if  B_element.get() == "" or B_key_name.get() == "" or B_key_energy.get() == "" :
        error("Certaines valeures obligatoires ne sont pas remplies")
    else:    
        if B_file_in.get() == "" :
            file_in = 'data_in.csv'
        else :
            file_in = B_file_in.get()
            
        if B_file_out.get() == "" :
            file_out = 'data_out'
        else :
            file_out = B_file_out.get()

        if file_in not in os.listdir():
            error("Le fichier d'entrée : ",file_in ,"n'a pas été trouvé")
        
        
        else:
            element = B_element.get().split(',')
            key_energy = B_key_energy.get()
            key_name = B_key_name.get()
            
            # Erreur sur la valeurs des paramètres
            key_file = np.loadtxt(file_in,dtype=str,delimiter=',')[0]
            list_error = []
            for key in [key_name,key_energy]:
                if key not in key_file :
                    list_error.append(key)
            
            if len(list_error) != 0 :
                txt = "Les colonnes suivantes n'ont pas été trouvées dans le fichier d'entrée >>"
                for i in list_error :
                    txt+=' '+i
                error(txt)
            else :
                for key in element:
                    if key not in key_file :
                        list_error.append(key)
            
                if len(list_error) != 0 :
                    txt = "Les éléments suivants ne sont pas dans la première ligne du fichier d'entrée  >>"
                    for i in list_error :
                        txt+=' "'+i+'"'
                    error(txt," ","Attention : les éléments doivent être séparé pas des virgules, sans espaces supplémentaites","pas de virgule à la fin")
                else :
                    flag = True
                    root.destroy() # Efface la valeur de Entry


#Interface graphique
flag = False
root = tk.Tk()
#Police

couleur = {
'gris' : '#bababa',
'noir' : '#000000'
}

#couleur = ['#9FF781','#fd9201','#da1512']
font_1 = ("arial",20)
font_2 = ("arial",15)
font_3 = ("arial",10)
# Création page
root.geometry('550x350')
root.config()
root.title("convexhull")


#Titre
title = tk.Label(root, text="Calcul d'enveloppe convexe", font=font_1)
title.pack()

# Lecture données d'entré


lf = tk.LabelFrame(root, text="Fichiers")
lf.place(x=5, y=40, height=70, width=540)

yrel = 0.25
#file_in
lentry = tk.Label(lf, text="Entrée (complet) :", font=font_3)
lentry.place(relx=0.02, rely=yrel)
B_file_in = tk.Entry(lf, font=font_3, width=17)
B_file_in.insert(0,'data_in.csv')
B_file_in.place(relx=0.23, rely=yrel)

#file_out
lentry = tk.Label(lf, text="Sortie (sans extention) :", font=font_3)
lentry.place(relx=0.48, rely=yrel)
B_file_out = tk.Entry(lf, font=font_3, width=17)
B_file_out.insert(0,'data_out')
B_file_out.place(relx=0.75, rely=yrel)



lf = tk.LabelFrame(root, text="Nom des colonnes dans le tableaux d'entrée")
lf.place(x=5, y=120, height=140, width=540)

#key_name
yrel = 0.1

lentry = tk.Label(lf, text="Nom* :", font=font_3)
lentry.place(relx=0.085, rely=yrel)
B_key_name = tk.Entry(lf, font=font_3, width=17)
B_key_name.place(relx=0.17, rely=yrel)

# key_energy
lentry = tk.Label(lf, text="Energie* :", font=font_3)
lentry.place(relx=0.565, rely=yrel)
B_key_energy = tk.Entry(lf, font=font_3, width=17)
B_key_energy.place(relx=0.68, rely=yrel)

yrel = 0.5

# element
lentry = tk.Label(lf, text="Références* (éléments purs) :", font=font_3)
lentry.place(relx=0.015, rely=yrel)
B_element = tk.Entry(lf, font=font_3, width=48)
xrel = 0.35
B_element.place(relx=xrel, rely=yrel)

lentry = tk.Label(lf, text="Chaque référence doit être séparé par une virgule", font=font_3)
lentry.place(relx=xrel, rely=0.7)

#Légedende
lentry = tk.Label(root, text="* Obligatoire", font=font_3)
lentry.place(relx=0.01, rely=0.75)

# Bouton
run = tk.Button(root, font=font_2, text="Run", command=clic_bouton) # Création du bouton
run.place(relx=0.45, rely=0.83)

# Affichage Fenêtre
root.mainloop()

#-----------------------------------------------------------------------------------------------------------------------------------------------------
if not flag:
    exit()

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