# convexhull
Le code `convexhull` permet de tracer l'enveloppe convexe à N dimension à partir de l'énergie d'un set de phases. Il est basé sur l'algotihme quickhull [1].
On peut obtenir la liste des phases sur l'enveloppe convexe, l'écart en énergie par rapport à l'enveloppe. On peut également obtenir la décomposition de phases à une composition donnée ou à partir d'une phases qui n'est pas sur l'enveloppe convexe.  
## Format
Le code est rédigé dans le fichier `cvh.py` sous forme d'une classe python `cvhull` avec plusieurs atributs et méthodes. La classe est commenté, on peut donc obtenir sont utilisation et celles de ses méthodes en appliquant la fonction help. 
Certaines fonctions sont répertoriés dans le fichier `cvh_func.py`.
Dans le dossier il y a un notebook jupyter d'exemple qui montre les utilisations principales du code (`exemple_cvh.pynb`)
## Paquets nécessaire à son utilisation
* cvh.py et cvh_func.py
* numpy
* scipy
## Aide de la classe `cvhull`
La classe prend en entrée un nombre N de phases de c éléments et caculs les phases stables et instables
### Paramètres :
* `nom` str[N] : liste des noms des phases d'entrée
* `composition` float[N,c] : liste des compsitions des phases
* `energie` float[N] : liste des énergie des phases
* `prec` int : nombre de chiffres après la virgule utlisé pour les compositions (défaut=4)
### Attributs :
* `cvhull.nom` str[N] : Liste des noms des phases d'entrée
* `cvhull.composition` np.array([N,c],float) : Liste des compsitions des phases
* `cvhull.energie` np.array([N],float)   : Liste des énergie des phases
* `cvhull.points` np.array([N,c],float) : Liste de points dans le format quickhull [composition indépendante + energie]

* `cvhull.npoints` int : Nombre de phases
* `cvhull.dim` int : Dimension des phases
* `cvhull.prec` int : Nombre de chiffres après la vigures considéré pour les compostions

* `cvhull.stable` np.array([N],int) : Liste des indexs des phases stables
* `cvhull.instable` np.array([N],int) : Liste des index des phases instable

* `cvhull.equations` np.array([n,c+1],float) : Liste des paramètres $p_i$ des équations des hyper-plans de l'enveloppe convexe (ex à 2D: $p_1x+p_2y+p_3 = 0$ )
* `cvhull.sommets` np.array([n,c+1],float) : Liste des sommets de l'hyper-plan
### Méthodes :
* `cvhull.energie_hull(compo)`\
Méthode qui donne l'énergie de l'enveloppe convexe à la composition `compo` (float[c]) donnée
* `cvhull.energie_compose(nom)`\
Méthode qui donne l'énergie de l'énergie de la phase `nom` (str)
* `cvhull.phases_stables(compo,prec=4)`\
Méthode qui donne la décomposition de phases à la composition `compo` (float[c]) donnée. les valeurs sont exprimées avec `prec` chifres apres la virgules (defaut :4)
* `cvhull.distance_hull(nom,prec=4,sortie=True)`\
Méthode qui donne l'écart en énergie par rapport à l'enveloppe convexe de phases `nom` (str).
Si elle n'est pas stable elle donne la décomposition de phases (peut être bloqué avec sortie=False).
Les valeur sont donnée avec `prec` chiffres après la virgule (defaut=4)

## Références
[1] C. B. Barber,D. P. Dobkin and H. Huhdanpaa, ACM Transactions on Mathematical Software, 22(4),469 (1996)
